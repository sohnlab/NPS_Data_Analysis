function [output_table, filtered_data, y_baseline, t_filtered, fs_filtered] = ...
    mNPS_procData(filepath, ch_height, De_np, wC, thresholds, sampleRate, ASLS_param, eventlength_filt)
    %   ========================================================
    %   2026/09/01 Chang
    %   Modified based on the published code in the JoVE paper
    %   ========================================================
    %
    %   Reads all mNPS data, analyzes data, and returns final output matrix.
    %   Needs a vector of 2 thresholds for initial thresholding. Afterwards,
    %       QC, thresholding, and analysis for individual pulses will need to
    %       proceed with user input.
    % 
    %   Inputs, outputs, and the output table columns are documented in README_v2.md.

    %% parse inputs

    if nargin < 5 || isempty(thresholds)
        thresholds = [1e-4, 1e-3];
        fprintf('Auto thresholds set to %3.2e, %3.2e\n',thresholds);
    end

    if nargin < 6 || isempty(sampleRate)
        sampleRate = 50000;
        fprintf('default sample rate used: %d Hz\n', sampleRate);
    end

    % default ASLS parameters
    if nargin<7 || isempty(ASLS_param)
        ASLS_param = struct();
        ASLS_param.lambda = 1e9;            % default 1e5; larger=smoother, smaller=wiggly-er (may not be unit-independent)
        ASLS_param.p = 3e-3;                % default 0.01; 0>p>1 (as low as possible while still converging)
        ASLS_param.noise_margin = 1e-4;     % default 0; allows baseline to sit within the baseline noise
        ASLS_param.max_iter = 20;           % default 5; just make sure it converges
    end

    if nargin<8 || isempty(eventlength_filt)
        eventlength_filt = 2000;
    end

    %% read all, search for pulses, get column information

    % load data
    load(filepath,'data');
    if ~( size(data,1)==1 && size(data,2)>2 )
        error('the variable `data` in filepath must be a 1xn array of doubles');
    end
    
    % First call readData: Scan whole file
    [all_out, ~, ~, outcols, outunits, rec_des, filtered_data, y_baseline, t_filtered, fs_filtered] = ...
        mNPS_readData(data, sampleRate, ch_height, De_np, wC, thresholds, false, false, ASLS_param);

    %% remove duplicate files to read
    [uni_win, output_matrix] = mNPS_cleanKim(all_out);

    clear all_out

    %% analyze all time windows, prompt user for input to see if the pulse looks good
    good_index = 1; % index through output of good pulses
    i = 1;
    new_th = thresholds; % reset threshold values
    
    % quieter command window while stepping through pulses:
    %   backtrace off drops the "> In mNPS_readData (line 41) ..." stack printed
    %   under every warning; the message itself still appears.
    % onCleanup puts both back when this function exits, including on error or
    % Ctrl-C, so the session is left as it was found.
    bt_state = warning('off', 'backtrace');
    cf_state = warning('off', 'curvefit:cfit:subsasgn:coeffsClearingConfBounds'); % annoying warning when fit fails to converge
    restore_warnings = onCleanup(@() warning([bt_state, cf_state]));

    % window highlighted in figure 42 and saved by [ 1 ]; [ 2 ] moves it.
    % sel_pulse records which pulse the selection belongs to, so it never
    % carries over to the next one.
    sel_window = 1;
    sel_pulse  = 0;

    while (i <= size(uni_win,1)) % <= : the last candidate is a pulse too

        % checking for working pulses
        searchflag  = true; % true to search, false to stop
        retryflag   = true; % true for first attempt, false if trying again with default thresholds

        while (searchflag)
            try
                fprintf('\nNow reading index: %d\n',uni_win(i));
                fprintf('Progress: %2.1f %%\n',i/length(uni_win)*100);

                % set window size parameters
                min_windowsize_filt = ceil(sampleRate/20); % for LPF padding, ensure window size cover >1sec
                startoffset_filt    = 200; % # filtered samples included before the first detected peak
                windowsize_filt     = eventlength_filt * 1.2; % # filtered samples w/ 20% buffer
                
                % set window indices
                startix_filt = uni_win(i)-startoffset_filt;
                endix_filt   = startix_filt + max(min_windowsize_filt,windowsize_filt);

                startix      = max(1, 20*startix_filt);             % ensure valid index
                endix        = min(length(data), 20*endix_filt);    % ensure valid index

                % extract iterdata window & check length
                iterdata = data(startix:endix);

                if length(iterdata) < sampleRate+1
                    warning('iterdata has only %u samples, but sampleRate=%u; this can cause problems with pad_data in LPF', length(iterdata), sampleRate);
                end
                
                % measure a new pulse
                % Second call readData: Threshold the pulses
                [iter_out, emptyflag, ~,~,~,~,~,~,~,~] = mNPS_readData(iterdata, sampleRate, ch_height, De_np, wC, new_th, false, false, ASLS_param);

                % iter_out:     output of one iteration
                % emptyflag:    skip pulse if TRUE
                % auto:         values for computing auto-threshold value

                % function arguments: read iterdata from uni_win entries
                % plot and fit off for speed
                if emptyflag
                    fprintf('No pulse! Skipping...\n');
                    i = i+1;
                    searchflag = true;
                else
                    % Third call readData: Get window range and plot
                    [iter_out, ~, auto, ~,~,~,~,~,~,~, win_range] = mNPS_readData(iterdata, sampleRate, ch_height, De_np, wC, new_th, true, true, ASLS_param);
                    searchflag = false;

                    % mNPS_readData always highlights window 1. Put the
                    % highlight back where [ 2 ] left it, as long as this is
                    % still the same pulse and that window survived the re-read.
                    if sel_pulse ~= i
                        sel_pulse  = i;
                        sel_window = 1;
                    elseif sel_window > size(iter_out,1) || iter_out(sel_window,1) == 0
                        sel_window = 1;
                    elseif sel_window ~= 1
                        mNPS_showWindow(sel_window);
                    end
                end

            catch ME

                %% print errors to command line for debug purposes
                fprintf('-----\n%s\n',ME.identifier);
                for errorstack_i = 1:length(ME.stack)
                    fprintf('Line: %d --- %s\n',ME.stack(errorstack_i).line,ME.stack(errorstack_i).name);
                end
                if retryflag % retry once with default thresholds
                    new_th = thresholds;
                    searchflag = true;
                    fprintf('Error occured, retrying with default thresholds...\n');
                    retryflag = false;
                else % fail on second try -> skip entirely
                    new_th = thresholds;
                    i = i + 1; % skip this pulse
                    fprintf('Error occured, skipping this file...\n');
                    retryflag = true;
                end
            end

            if i > size(uni_win,1) % return if finished
                break
            end
        end

        % The search above can run i past the last window (empty/error pulses
        % increment i). That only broke the inner loop, so guard here too:
        % otherwise the prompt shows a stale pulse and saving indexes
        % uni_win(i) out of range.
        if i > size(uni_win,1)
            break
        end

        % prompt for input to determine next operation
        fprintf(['\nIs anything wrong?\n' ...
                    '  ENTER  skip this pulse\n' ...
                    '  [ 1 ]  Save window %d (highlighted)\n' ...
                    '  [ 2 ]  Adjust window\n' ...
                    '  [ 3 ]  Auto threshold\n' ...
                    '  [ 4 ]  Adjust thresholds\n' ...
                    '  [ X ]  Stop and Output\n' ...
                    ], sel_window);
        fprintf('\n%d cell(s) saved\n',good_index-1);
        fprintf('Current thresholds: %3.2e, %3.2e\n',new_th);
        OK = strtrim(input('---\n','s'));

        %% input case structures
        %   the legacy keys are kept as undocumented aliases so that existing muscle memory still works
        switch OK
            
            % empty input: skip the pulse
            case ''

            fprintf('Skipping this pulse...\n');
            i = i+1;
            new_th = thresholds;

            % save the highlighted window
            case {'1', 'P','p', '.', '/'}

            % make and display table of pulses and indices.
            %   Numbered the same way as the '2' branch, the plot labels and
            %   win_range: one entry per row of iter_out. The old version
            %   filtered on start index > 100 and renumbered what was left,
            %   so its 'Window 1' could be a different window from the one
            %   highlighted in figure 42 and saved below.
            indices = iter_out(:,1);
            winds = 1:length(indices);
            table_data = [winds', indices];

            % clean up empty table entries
            cci = 1;
            stopc = size(table_data,1);
            while(cci <= stopc)
                if table_data(cci,2) == 0
                    table_data(cci,:) = [];
                    stopc = stopc - 1;
                else
                    cci = cci + 1;
                end
            end

            if ~isempty(table_data) % make sure WindowTable is NOT empty

                disp( array2table(table_data,'VariableNames',{'Window','StartIndex'}) );
                iter_out_index = sel_window; % the window highlighted in figure 42
                fprintf('Ok, saving data...\n');

                output_matrix(good_index,:) = iter_out(iter_out_index,:); % save to output matrix
                output_matrix(good_index,1) = output_matrix(good_index,1) + uni_win(i) - 200;
                good_index = good_index + 1;
                i = i+1;
                new_th = thresholds; % reset thresholds

            else
                fprintf('Pulse table is empty; force retry\n');
                new_th = thresholds; % reset thresholds
            end

            % user picks window
            case {'2', '//'}

            % make and display table of pulses and indices
            indices = iter_out(:,1);
            winds = 1:length(indices);
            table_data = [winds', indices];

            % clean up empty table entries
            cci = 1;
            stopc = size(table_data,1);
            while(cci <= stopc)
                if table_data(cci,2) == 0
                    table_data(cci,:) = [];
                    stopc = stopc - 1;
                else
                    cci = cci + 1;
                end
            end

            if ~isempty(table_data) % make sure WindowTable is NOT empty

                fprintf('Window %d is highlighted (the one that [ 1 ] would save).\n', sel_window);

                while true

                    disp( array2table(table_data,'VariableNames',{'Window','StartIndex'}) );
                    iter_out_index_str = strtrim(input([ ...
                        'Select window to highlight:\n' ...
                        '  ENTER  back to the main menu\n' ...
                        '  [ # ]  window number from the table above\n' ...
                        '---\n'],'s'));

                    % leave selection mode with the highlight where it is
                    if isempty(iter_out_index_str) || strcmp(iter_out_index_str,'0')
                        fprintf('Leaving window selection, back to this pulse.\n');
                        break
                    end

                    iter_out_index = str2double(iter_out_index_str);
                    if isnan(iter_out_index) || ~isscalar(iter_out_index) || ...
                            ~ismember(round(iter_out_index), table_data(:,1))
                        fprintf('Unrecognized input.\n');
                        beep
                        continue
                    end
                    iter_out_index = round(iter_out_index);

                    % move the highlight onto this window, re-coloring the
                    % segments and redrawing the fit panel
                    sel_window = iter_out_index;
                    mNPS_showWindow(sel_window);
                    fprintf(['Window %d spans index %d-%d (shaded).\n' ...
                             'red = sizing, yellow = contraction, blue = recovery.\n' ...
                             '[ 1 ] on the main menu now saves this window.\n'], ...
                        sel_window, win_range(sel_window,1), win_range(sel_window,2));
                    break

                end

            else
                fprintf('Pulse table is empty; force retry\n');
                new_th = thresholds; % reset thresholds
            end

            % auto threshold
            case {'3', '+','=', ''''}

            new_th = auto_both(auto, thresholds);
            sel_window = 1; % new thresholds renumber the windows
            fprintf('Auto-set both: bottom %3.2e, top %3.2e\n', new_th(1), new_th(2));

            % threshold-adjustment mode. 
            case {'4', 'T','t', 'Y','y'}
            
            % local call to update the new threshold
            reread = @(th) mNPS_readData(iterdata, sampleRate, ch_height, De_np, wC, th, true, true, ASLS_param);
            [new_th, thresholds, iter_out, auto, win_range] = ...
                threshold_menu(new_th, thresholds, reread, iter_out, auto, win_range);
            sel_window = 1; % new thresholds renumber the windows

            % unrecognized input: stop analyzing data and save processed events so far
            otherwise

            if confirm_stop_request('Unrecognized input.')
                break
            end
        end
    end

    %% finish up

    % remove empty rows from the output arrray
    output_data  = output_matrix( ~all(output_matrix==0,2), :);
    
    % don't include recovery time twice
    output_data  = output_data(:, ~strcmp(outcols, 'T_rec'));
    column_names = outcols(~strcmp(outcols, 'T_rec'));
    column_units = outunits(~strcmp(outcols, 'T_rec'));

    % description for the recovery category
    column_descriptions = cell(size(column_names));
    column_descriptions(:) = {''};
    column_descriptions(strcmp(column_names, 'rec_cat')) = {rec_des};
    
    % convert output array to a table with variable names & info
    output_table = array2table(output_data, 'variablenames', column_names);
    output_table.Properties.DimensionNames{2} = 'cell_data';
    output_table.Properties.VariableUnits = column_units;
    output_table.Properties.VariableDescriptions = column_descriptions;

    % play an alert sound to signal that processing is finished
    % (guarded: a missing/busy audio device must not abort processing before
    % the output table is finalized and returned below)
    try
        t = linspace(0,1,2^16); % time-samples
        Fs = 2^16; % sampling frequency
        y = 0.3*exp(-4*t).*sin(t*2*pi*440); % sinusoid
        sound(y,Fs);
    catch
        % audio device unavailable (headless/remote/busy) - skip the beep
    end

    % remove duplicate cell detections (keep last-recorded instance)
    output_table = remove_duplicate_rows(output_table);

    % second wCDI calculation, normalized by population mean Uflow.
    % wCDI_original (from mNPS_readData) divides each cell by its own V_np
    U_flow = mean(output_table.V_np, 'omitnan'); % [mm/s]
    output_table.wCDI_mean_Uflow = ...
        output_table.V_c ./ U_flow .* output_table.diameter ./ ch_height;

    % label both wCDI columns, so the table itself says which normalization is
    % which. (Adding a variable pads VariableUnits/VariableDescriptions with ''.)
    orig_col = strcmp(output_table.Properties.VariableNames, 'wCDI_original');
    new_col  = strcmp(output_table.Properties.VariableNames, 'wCDI_mean_Uflow');
    output_table.Properties.VariableUnits{new_col} = 'dimensionless';
    output_table.Properties.VariableDescriptions{orig_col} = ...
        'wCDI normalized by this cell''s own V_np (JoVE convention)';
    output_table.Properties.VariableDescriptions{new_col} = ...
        'wCDI normalized by U_flow, the mean V_np of every accepted cell in this run';

    fprintf('U_flow = %.3f mm/s over %d cell(s)\n', U_flow, height(output_table));

    fprintf('Done reading, check output!\n');

end

function [new_th, thresholds, iter_out, auto, win_range] = ...
        threshold_menu(new_th, thresholds, reread, iter_out, auto, win_range)
    %THRESHOLD_MENU dedicated mode for adjusting the two detection thresholds.
    %   Every change is applied immediately: the pulse is re-read and figure 42
    %   redrawn, without going back to the main menu.
    %
    %   new_th     = [bottom, top] currently in use
    %   thresholds = the defaults new_th is reset to at every new pulse, and
    %                whose bottom value floors the auto-set bottom; [ 6 ]
    %                replaces them, so the change carries to later pulses
    %   reread     = handle taking a threshold pair and returning a fresh
    %                mNPS_readData result for the pulse on screen
    %   iter_out, auto, win_range = the current results, replaced on each re-read

    while true

        fprintf(['\nAdjust thresholds  (bottom %3.2e, top %3.2e)\n' ...
                 '                 (default bottom %3.2e, top %3.2e)\n' ...
                 '  [ 1 ]  auto-set both\n' ...
                 '  [ 2 ]  auto-set top only\n' ...
                 '  [ 3 ]  auto-set bottom only\n' ...
                 '  [ 4 ]  enter top manually\n' ...
                 '  [ 5 ]  enter bottom manually\n' ...
                 '  [ 6 ]  make the current values the default\n' ...
                 '  ENTER  back to the main menu\n'], ...
                 new_th(1), new_th(2), thresholds(1), thresholds(2));
        choice = strtrim(input('---\n','s'));

        switch choice

            case ''
                return

            case '1'
                new_th = auto_both(auto, thresholds);
                fprintf('Auto-set both: bottom %3.2e, top %3.2e\n', new_th(1), new_th(2));

            case '2'
                new_th(2) = auto*0.85;
                fprintf('Auto-set top: %3.2e\n', new_th(2));

            case '3'
                new_th(1) = max(new_th(2)*0.12, thresholds(1));
                fprintf('Auto-set bottom: %3.2e\n', new_th(1));

            case '4'
                th_input = prompt_threshold('top');
                if isnan(th_input), continue, end
                new_th(2) = th_input;
                fprintf('New top threshold: %3.2e\n', new_th(2));

            case '5'
                th_input = prompt_threshold('bottom');
                if isnan(th_input), continue, end
                new_th(1) = th_input;
                fprintf('New bottom threshold: %3.2e\n', new_th(1));

            case '6'
                % the defaults are what every later pulse starts from, so this
                % is the only change here that outlives the pulse on screen.
                thresholds = new_th;
                fprintf('New default thresholds: bottom %3.2e, top %3.2e\n', ...
                    thresholds(1), thresholds(2));
                continue % nothing changed for this pulse; no need to re-read

            otherwise
                fprintf('Unrecognized input.\n');
                beep
                continue

        end

        if new_th(1) >= new_th(2)
            fprintf('WARNING: the bottom threshold is not below the top one.\n');
        end

        % apply straight away so the effect is visible without leaving here.
        % A threshold that finds nothing (or errors) leaves the last good
        % result in place, so the auto values stay usable.
        try
            [new_out, emptyflag, new_auto, ~,~,~,~,~,~,~, new_range] = reread(new_th);
            if emptyflag
                fprintf('No pulse detected with these thresholds.\n');
            else
                iter_out = new_out; auto = new_auto; win_range = new_range;
                mNPS_showWindow(1);
            end
        catch ME
            fprintf('Could not re-read with these thresholds: %s\n', ME.message);
        end

    end

end

function new_th = auto_both(auto, thresholds)
%AUTO_BOTH top threshold from the auto value, bottom scaled off it and
%   floored at the original default.

    new_th(2) = auto*0.85;
    new_th(1) = max(new_th(2)*0.12, thresholds(1));

end

function th = prompt_threshold(which_th)
%PROMPT_THRESHOLD read one threshold from the user; NaN means "no change".
%   str2double returns NaN for anything that isn't a number, so the value is
%   checked here: a NaN or non-positive threshold silently disables all
%   thresholding (abs(x) < NaN is never true) and every sample then reads as
%   a pulse.

    th_str = strtrim(input(sprintf([ ...
        'New %s threshold:\n' ...
        '  ENTER  cancel, leave it unchanged\n' ...
        '  [ # ]  a positive number, e.g. 5e-4\n' ...
        '---\n'], which_th), 's'));

    if isempty(th_str)
        th = NaN;
        return
    end

    th = str2double(th_str);
    if ~isfinite(th) || th <= 0
        fprintf('Not a valid threshold, no change made.\n');
        beep
        th = NaN;
    end

end

function should_stop = confirm_stop_request(reason)
    if nargin < 1 || isempty(reason)
        reason = 'Unrecognized input.';
    end

    fprintf('%s\n', reason);
    confirm_stop = input([ ...
        'Stop analyzing and output what is saved?\n' ...
        '  ENTER  no, keep going\n' ...
        '  [ Y ]  yes, stop and output\n' ...
        '---\n'], 's');
    confirm_stop = strtrim(confirm_stop);
    if isempty(confirm_stop)
        should_stop = false;
    else
        should_stop = any(strcmpi(confirm_stop, {'y','yes'}));
    end

    if should_stop
        fprintf('Stopping and dumping data.\n');
    else
        fprintf('Continuing without stopping.\n');
    end
end

%% remove obvious duplicates (keep last-recorded instance)

function new_table = remove_duplicate_rows(old_table)
    [~,ia,~] = unique(old_table.start_ix, 'last');
    new_table = old_table(ia,:);
end


