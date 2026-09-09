function [OUT_array, empty, auto_thresh_value, column_names, column_units, rec_cat_description, filtered_data, y_baseline, t_filtered, fs_filtered, win_range] = ...
    mNPS_readData(data_vector, sampleRate, ch_height, De_np, wC, thresholds, plotflag, fitflag, ASLS_param)
    %   ========================================================
    %   2026/09/01 Chang
    %   Modified based on the published code in the JoVE paper
    %   ========================================================
    %
    % Reads mNPS data and returns OUT_array matrix
    % 
    % Inputs, outputs, and the output table columns are documented in README_v2.md.

    if nargin<9
        ASLS_param = []; % use defaults
    end

    %% SECTION 0: device parameters

    % Segment layout
    total_segs   = 6;
    num_ref_segs = 2;
    num_rec_segs = 3;
    % check layout
    if total_segs ~= (num_ref_segs + num_rec_segs + 1)
        error('check device geometry! number of segments doesn''t match!');
    end

    % Channel Geometry [um]
    L       = 6750;             % total length of the NPS channel (from start of 1st pore to end of last pore)
    npL_ref = [700, 700];       % sizing pore length array
    wNP     = 22;               % sizing pore width
    sqL     = 3000;             % contraction pore length
    npL_rec = [700, 700, 700];  % recovery pore length array

    % calculate De_c based on De_np
    De_c = De_np*(wC/wNP)^(0.5); % D_e (effective diameter) for contraction segment

    %% SECTION 1: load data and perform basic signal conditioning

    Fs = sampleRate/1000;       % convert to kHz

    N  = 20;                    % downsample factor

    % flip if negative
    if (sum(data_vector<0) / length(data_vector)) > 0.5
        data_vector = -data_vector;
        warning('inverting raw data because the majority of datapoints were negative');
    end

    y_smoothed = fastsmooth(data_vector',200,1,1);  % perform rectangular smoothing

    ym = downsample(y_smoothed,N);                  % downsample by period N

    if size(ym,1) < size(ym,2)
        ym = ym';                                   % transpose if vector is of the wrong dimension
    end

    % low-pass filter
    %   - default bw_lpf would be 5/min_pulse_width [Hz]
    %   - observed short pulse in 0220617_A549_dev2B_w12_p25_try1.mat was 133 samples at (50kHz / 20) >>> 0.0532sec
    bw_lpf        = 50;
    fs_filtered   = sampleRate/N;
    filtered_data = unpad_data(lowpass(pad_data(ym-ym(1), fs_filtered), bw_lpf, fs_filtered), fs_filtered) + ym(1);

    % remove baseline
    ym = filtered_data;
    y_baseline = -1 * ASLS(-1*ym, ASLS_param);
    y_detrend = ym - y_baseline;

    % generate downsampled time vector
    t_raw = (0:length(data_vector)-1) / sampleRate;
    t_filtered = (0:length(filtered_data)-1).' / fs_filtered + t_raw(1); % [sec]

    %% SECTION 2: threshold signal by differences
    % take the difference of ym, threshold by lower value thresholds for differences, user provided

    ym_diff = diff(y_detrend);                  % compute difference

    ym_diff(abs(ym_diff) < thresholds(1)) = 0;  % threshold values below thresholds(1)
    ym_diff(1:10) = 0;                          % zero out the first few values
    ym_diff(end-10:end) = 0;                    % zero out the last few values

    % ensure all values in squeeze channel are zero.
    %   Only the first squeeze in the vector is blanked, matching the JoVE
    %   behaviour. On a whole-recording call that means later contractions stay
    %   in ym_diff; the per-event calls see one contraction each, so for them
    %   this is the whole story.
    squeeze_begin = find(ym_diff <= -thresholds(2), 1); % find where squeeze starts

    if ~isempty(squeeze_begin)

        % advance by 1 until large difference ends. Bounded: a squeeze running
        % to the end of the vector must not index past it.
        while squeeze_begin < numel(ym_diff) && ym_diff(squeeze_begin) <= -thresholds(2)
            squeeze_begin = squeeze_begin + 1;
        end

        squeeze_end = find(ym_diff(squeeze_begin+20:end) >= thresholds(2), 1); % find end of squeeze
        if ~isempty(squeeze_end)
            % zero out the squeeze channel
            zero_last = min(squeeze_begin + squeeze_end + 20, numel(ym_diff));
            ym_diff(squeeze_begin+1 : zero_last) = 0;
        end

    end

    %% SECTION 3: identify nonzero differences (nz_mat is the matrix of nonzero differences)

    nz_mat = ones(2,length(nonzeros(ym_diff))); % preallocation
    k = 1; % index for A
    for i = 1:length(ym_diff)
        if (ym_diff(i) ~= 0) % look for nonzeros
            nz_mat(1,k) = i; % array index of nonzero
            nz_mat(2,k) = ym_diff(i); % nonzero value
            k = k+1;
        end
    end

    %% SECTION 4: remove error from A

    k=1;
    while (k < ceil(log(length(nz_mat))))
        i=1;
        while (i < length(nz_mat))

            % Case 1: current and next both positive && next > current
            if nz_mat(2,i) > 0 && nz_mat(2,i+1) > 0 && ...
                    nz_mat(2,i+1) > nz_mat(2,i)
                % move next into current
                nz_mat(:,i) = nz_mat(:,i+1);

            % Case 2: current and next both positive && current > next
            elseif nz_mat(2,i) > 0 && nz_mat(2,i+1) > 0 && ...
                    nz_mat(2,i) > nz_mat(2,i+1)
                % move current into next
                nz_mat(:,i+1) = nz_mat(:,i);

            % Case 3: current and next both negative && current > next
            elseif nz_mat(2,i) < 0 && nz_mat(2,i+1) < 0 && ...
                    nz_mat(2,i) > nz_mat(2,i+1)
                % move next into current
                nz_mat(:,i) = nz_mat(:,i+1);

            % Case 4: current and next both negative && next > current
            elseif nz_mat(2,i) < 0 && nz_mat(2,i+1) < 0 && ...
                    nz_mat(2,i+1) > nz_mat(2,i)
                % move current into next
                nz_mat(:,i+1) = nz_mat(:,i);
            end

            i=i+1;
        end
        k=k+1;
    end

    %% SECTION 5: remove repeats in A

    unique_xs = unique(nz_mat(1,:)); % unique values
    unique_is = ones(1,length(unique_xs)); % unique indices

    for i = 1:length(unique_is)
        unique_is(i) = find(nz_mat(1,:) == unique_xs(i), 1);
    end
    unique_ys = nz_mat(2,unique_is);
    nz_mat = [unique_xs; unique_ys];

    %% SECTION 6: rectangularize pulses

    % This pass only ever reads y_detrend and nz_mat, and modifies neither, 
    % reduce the original 50 repeats to 1.
    ym_rect = y_detrend;
    i = 1;
    while (i < size(nz_mat,2))

        if nz_mat(2,i) < 0 && nz_mat(2,i+1) > 0 % look for sign change in differences
            % replace all values in between with mean
            ym_rect(nz_mat(1,i):nz_mat(1,i+1)) = ...
                mean(y_detrend(nz_mat(1,i):nz_mat(1,i+1)));
        end

        i=i+1;
    end

    %% SECTION 7: Plot figures if flag is true

    if max(ym_diff) < thresholds(1) % waste of time, no pulse
        plotflag = false;
        empty = true;
    else
        empty = false;
    end

    if plotflag

        if isempty( findobj('type','figure', 'number',42) ) % create a new figure
            Pix_SS = get(0,'screensize');
            figh = figure(42);
            figsize = [0.1 0.1 0.45 0.75]*Pix_SS(4);
            set(figh,'units','pixels','pos',figsize);
        else % use the existing figure (don't change size/location)
            figh = figure(42);
        end
        
        % take top and bottom 3 values
        nsorted_d = sort(ym_diff);
        min_vals = nsorted_d(1:3);

        psorted_d = sort(ym_diff,'descend');
        max_vals = psorted_d(1:3);

        % set auto-thresholds
        if abs(min_vals(3)) < abs(max_vals(3))
            auto_thresh_value = abs(min_vals(3));
        else
            auto_thresh_value = abs(max_vals(3));
        end

        % three signal panels sharing one x-axis, plus a recovery-fit panel.
        % Positions are set explicitly (rather than via subplot) to keep the
        % stack tight; they are re-applied after plotting because newplot
        % resets most axes properties.
        clf(figh);
        ax1 = axes(figh); ax2 = axes(figh); ax3 = axes(figh); ax4 = axes(figh);

        % difference plot
        plot(ax1, ym_diff,'k-', 'LineWidth',1);
        title(ax1,'y_{diff}');
        grid(ax1,'on');
        ax1.XMinorGrid = 'on';
        diff_pad = 0.05 * max(max(ym_diff)-min(ym_diff), eps);
        axis(ax1, [0, length(ym_diff), min(ym_diff)-diff_pad, max(ym_diff)+diff_pad]);

        % plot thresholds
        hold(ax1,'on');
        for yi=[1,-1]
            yline(ax1, yi*thresholds(1), 'b', 'linew',1);
            yline(ax1, yi*thresholds(2), 'b--', 'linew',1.5);
        end
        hold(ax1,'off');

        % rectangularized
        plot(ax2, ym_rect, 'k', 'linew',1);
        title(ax2,'y_{rect}');
        grid(ax2,'on');
        ax2.XMinorGrid = 'on';
        % pad the pulse depth by 5% instead of a fixed top limit, so shallow
        % pulses are not squashed against the bottom of the panel
        rect_pad = 0.05 * max(max(ym_rect)-min(ym_rect), eps);
        axis(ax2, [0, length(ym_rect), min(ym_rect)-rect_pad, max(ym_rect)+rect_pad]);
        xlabel(ax2, 'filtered sample index'); % bottom of the signal stack

        % smoothed
        plot(ax3, ym, 'k', 'linew',1);
        hold(ax3,'on'); plot(ax3, ym-y_detrend, 'm', 'linew',1); hold(ax3,'off');
        title(ax3,'y_{LP}');
        grid(ax3,'on');
        ax3.XMinorGrid = 'on';
        % fit both traces, with a strip of headroom on top for the window
        % number labels that Section 8 places just above the data
        y_base = ym - y_detrend;
        lp_lo = min([ym(:); y_base(:)]);
        lp_hi = max([ym(:); y_base(:)]);
        lp_range = max(lp_hi-lp_lo, eps);
        axis(ax3, [0, length(ym), lp_lo-0.05*lp_range, lp_hi+0.22*lp_range]);

        % recovery-fit panel: filled in by mNPS_showWindow once the fits exist
        title(ax4, 'recovery fit');

        % apply the tight stack and drop the duplicated x tick labels.
        % Visual order, top to bottom: y_diff, y_LP, y_rect, recovery fit.
        % Only the bottom signal panel keeps its tick labels and xlabel.
        L = 0.10; Wd = 0.86; h_sig = 0.20; h_fit = 0.15; gap = 0.042; b0 = 0.065;
        b3 = b0 + h_fit + gap + 0.02;
        b2 = b3 + h_sig + gap;
        b1 = b2 + h_sig + gap;
        set(ax1,'Position',[L b1 Wd h_sig], 'FontSize',9, 'XTickLabel',[]); % y_diff
        set(ax3,'Position',[L b2 Wd h_sig], 'FontSize',9, 'XTickLabel',[]); % y_LP
        set(ax2,'Position',[L b3 Wd h_sig], 'FontSize',9);                  % y_rect
        set(ax4,'Position',[L b0 Wd h_fit], 'FontSize',9);

        linkaxes([ax1,ax2,ax3], 'x');

    else
        auto_thresh_value = [];
    end

    %% SECTION 8: Detect NPS pulses
    % pulse_series is a matrix with the indices and parameters for rectangular pulses

    i=1;
    k = 0;
    backset = 10;
    pulse_series = ones(length(nz_mat),5);
    while (i < length(nz_mat))
        if nz_mat(2,i) < 0 && nz_mat(2,i+1) > 0 % starts negative and flips sign
            k = k + 1;
            pulse_series(k,1) = nz_mat(1,i); % Start index
            pulse_series(k,2) = nz_mat(1,i+1); % End index
            pulse_series(k,3) = mean(ym((nz_mat(1,i)-backset):(nz_mat(1,i)-backset+10))); % normalized baseline current
            pulse_series(k,4) = mean(y_detrend(nz_mat(1,i)+1:nz_mat(1,i+1)-1)); % avg current drop between pulses
            pulse_series(k,5) = std(y_detrend(nz_mat(1,i)+1:nz_mat(1,i+1)-1)); % std dev of current drop
        end
        i=i+1;
    end

    % remove empty entries in P
    cci = 1;
    stopc = size(pulse_series,1);
    while(cci <= stopc)
        if (pulse_series(cci,1) == 1 && pulse_series(cci,2) == 1)
            pulse_series(cci,:) = [];
            stopc = stopc - 1;
        else
            cci = cci + 1;
        end
    end

    % index range spanned by each candidate window: window k covers segments k
    % through k+total_segs-1, so it runs from the start of segment k to the end
    % of segment k+total_segs-1. Rows line up 1:1 with the rows of OUT_array and
    % with the 'Window' column of the selection table in mNPS_procData.
    num_windows = size(pulse_series,1) + 1 - total_segs;
    if num_windows >= 1
        win_range = [pulse_series(1:num_windows,1), pulse_series(total_segs:end,2)];
    else
        num_windows = 0;
        win_range = zeros(0,2);
    end

    % plot detected pulses
    if plotflag
        nstarts = pulse_series(:,1);
        nstops = pulse_series(:,2);
        hold([ax1,ax2,ax3], 'all');
    
        % on the diff plot
        plot(ax1, nstarts,ym_diff(nstarts), 'bo', 'linew',1.5, 'markersi',8);
        plot(ax1, nstops,ym_diff(nstops), 'ro', 'linew',1.5, 'markersi',8);

        % draw each detected segment on the rect and data plots. The colors
        % are set by mNPS_showWindow according to the role each segment
        % plays in the highlighted window; the handles are stashed on the
        % figure so mNPS_procData can re-color when a different window is
        % picked. Window 1 is the default, matching the default highlight.
        seg_lines = gobjects(length(nstarts), 2);
        for ii=1:length(nstarts)
            seg_ix = nstarts(ii):nstops(ii);
            seg_lines(ii,1) = plot(ax2, seg_ix, ym_rect(seg_ix), 'linew',3);
            seg_lines(ii,2) = plot(ax3, seg_ix, ym(seg_ix), 'linew',2.5);
        end

        W.lines = seg_lines;
        W.num_ref = num_ref_segs;
        W.num_rec = num_rec_segs;
        W.win_range = win_range;
        W.ax_sig = [ax1, ax2, ax3];
        W.ax_fit = ax4;
        W.rT = []; W.rdI = []; W.p1 = []; W.p2 = []; % filled in after Section 9
        setappdata(figh, 'mNPS_window_data', W);

        % label each candidate window at its start index, so the numbers match
        % the 'Window' column of the '//' selection table in mNPS_procData
        if num_windows >= 1
            % two staggered rows just above the traces, so that windows with
            % nearby start indices don't overprint each other. Anchored to the
            % top of the data rather than to a fraction of the axis, so the
            % labels stay close to it whatever the pulse depth.
            ytext = lp_hi + [0.16, 0.07]*lp_range;
            for ii = 1:num_windows
                xline(ax3, win_range(ii,1), ':', 'color',[0.5 0.5 0.5], 'linew',1);
                text(ax3, win_range(ii,1), ytext(1+mod(ii,2)), sprintf(' %d',ii), ...
                    'FontSize',8, 'FontWeight','bold', 'Color',[0.25 0.25 0.25], ...
                    'HorizontalAlignment','left', 'VerticalAlignment','middle');
            end
        end

        hold([ax1,ax2,ax3], 'off');
    end


    %% SECTION 9: Extract mNPS pulse data

    % preallocate array of NPS pulse data
    out = nan(length(pulse_series) - (total_segs-1), 12);

    % preallocate array of reference segment dT_np_segs
    dT_np_segs = nan(size(out,1), num_ref_segs);

    % preallocate the per-window recovery curve data (for the fit panel)
    rT_all = nan(size(out,1), num_rec_segs);
    rdI_all = nan(size(out,1), num_rec_segs);

    % column metadata and the recovery description are constants. They used to
    % be assigned inside the loop below, so a vector with fewer than total_segs
    % detected segments left them undefined and the function failed with
    % "Output argument not assigned" instead of returning an empty result.
    out_cols = {'start_ix', 'I_baseline', 'dI_np', 'dI_c', 'dI_c_std', ...
        'dT_np', 'dT_c', 'T_rec', 'fo_p1', 'fo_p2', 'gof_rsquare', 'rec_cat'};
    out_units = {'index', 'data units', 'data units', 'data units', 'data units', ...
        'ms',    'ms',   'ms',    '',      '',      '',        'categorical'};
    rec_cat_description = '0 = instant, 1-2 = transient, 3 = prolonged';

    for k = 1 : length(pulse_series)+1-total_segs

        % get relevant row numbers
        sq_k = k + num_ref_segs; % contraction segment
        ref_k_start = k; % first reference segment
        ref_k_end = sq_k - 1; % last reference segment
        rec_k_start = sq_k + 1; % first recovery segment
        rec_k_end = rec_k_start + num_rec_segs - 1; % last recovery segment
        % check segment indices
        if rec_k_end ~= (k + total_segs - 1)
            error('check segment indexing!');
        end

        start_index = pulse_series(k,1); % starting index
        I_baseline = pulse_series(k,3); % baseline current
        
        % average dI & dT in reference segments

        % average node-pore current drop in reference segments
        dI_np = -mean(pulse_series(ref_k_start:ref_k_end,4));

        % average node-pore transit time in reference segments [ms]
        %   *** not applicable in JOVE device designs (sNPS_ver2.1) bc the segments are of unequal lengths
        % dT_np = nan;

        % node-pore transit time in each reference segment [ms] (row vector)
        dT_nps = ( pulse_series(ref_k_start:ref_k_end,2) - pulse_series(ref_k_start:ref_k_end,1) )' ./Fs.*N;
        dT_np = mean(dT_nps,2);

        % dI & dT in contraction (sqeeze) segment
        dI_c = -pulse_series(sq_k,4); % squeeze current drop
        dI_c_std = pulse_series(sq_k,5); % std. dev. of squeeze current drop
        dT_c = (pulse_series(sq_k,2) - pulse_series(sq_k,1)) /Fs*N; % squeeze transit time (ms)

        % post-squeeze NP current drops
        dI_rec = -pulse_series(rec_k_start:rec_k_end, 4);
        
        %% determine recovery time & category
        % recovery time is determined when post-squeeze NP current drop
        %   reaches pre-squeeze NP current drop (within 8% error threshold)
        % if the cell has "instant" recovery, recovery time is defined as
        %   elapsed time between the end of the squeeze segment and the
        %   beginning of the first recovery segment
        % if the cell has "transient" recovery, recovery time is defined as
        %   elapsed time between the end of the squeeze segment and the
        %   beginning of the first segment where the cell was recovered
        
        % "recovered" is when dI_rec comes within 8% of dI_np or higher
        rec_tol = 0.08;
        
        if num_rec_segs ~= 3
            warning('recovery time & category are hard-coded for devices with 3 recovery segments!');
        end
        
        % cell was already recovered by the first recovery segment
        if (dI_np-dI_rec(1))/dI_np < rec_tol
            T_rec = (pulse_series(rec_k_start,1) - pulse_series(sq_k,2)) /Fs*N;
            rec_cat = 0;

        % cell didn't recover until the second recovery segment
        elseif (dI_np-dI_rec(2))/dI_np < rec_tol
            T_rec = (pulse_series(rec_k_start+1,1) - pulse_series(sq_k,2)) /Fs*N;
            rec_cat = 1;

        % cell didn't recover until the third recovery segment
        elseif (dI_np-dI_rec(3))/dI_np < rec_tol
            T_rec = (pulse_series(rec_k_start+2,1) - pulse_series(sq_k,2)) /Fs*N;
            rec_cat = 2;

        % by the third recovery segment, cell still hadn't recovered
        else
            T_rec = Inf;
            rec_cat = 3;

        end

        %% perform mNPS-r recovery curve fitting
        if fitflag

            if num_rec_segs ~= 3
                error('mNPS-r recovery fitting is hard-coded for devices with exactly 3 recovery segments');
            end

            % time-point of each recovery pulse: the START index of that
            % recovery segment, referenced to the END of the contraction
            % segment. Every point therefore uses the same index of its own
            % region and the same origin - the moment deformation stopped -
            % which is also the convention T_rec uses above.
            rT = ( pulse_series(rec_k_start:rec_k_end, 1) - pulse_series(sq_k,2) )' ...
                ./ Fs .* N; % in ms

            % populate vector rdI with current drop-amplitudes of recovery
            % pulses. Anticipate that rdI should increase
            %
            % NOTE: cannot approximate size for ellipsoid particle
            rdI = dI_rec';

            % use MATLAB fit() to fit a linear polynomial to recovery data.
            % fo has fields pertaining to fit parameters:
            %   fo.p1 is the slope of the line, not sign-bounded
            %   fo.p2 is the y-intercept, not sign-bounded
            % gof has fields pertaining to goodness of fit, but only
            % gof.rsquare will be used
            [fo, gof] = fit(rT',rdI','poly1');

        else
            fo.p1 = 0;
            fo.p2 = 0;
            gof.rsquare = 0;
        end

        out(k,:) = [start_index, I_baseline, dI_np, dI_c, dI_c_std, dT_np, dT_c, T_rec, fo.p1, fo.p2, gof.rsquare, rec_cat];

        % save transit time in each reference segment
        dT_np_segs(k,:) = dT_nps;

        % save the recovery curve for this window
        if fitflag
            rT_all(k,:) = rT;
            rdI_all(k,:) = rdI;
        end

    end

    % hand the per-window recovery curves to the fit panel of figure 42 and
    % show window 1, the default that mNPS_procData would save
    if plotflag
        W = getappdata(figh, 'mNPS_window_data');
        W.rT = rT_all;
        W.rdI = rdI_all;
        W.p1 = out(:,9);
        W.p2 = out(:,10);
        setappdata(figh, 'mNPS_window_data', W);
        mNPS_showWindow(1);
    end

    %% Section 9b: compute derived values

    calculated = zeros(size(out,1),9);

    % diameter (um)
    calculated(:,1) = ((out(:,3)./out(:,2)*De_np^2*L)./ ...
        (1+0.8*L/De_np*out(:,3)./out(:,2))).^(1/3);

    % strain (dimensionless)
    calculated(:,2) = (calculated(:,1)-wC)./calculated(:,1);

    % np velocity (mm/s = µm/ms)
    calculated(:,3) = mean(npL_ref./dT_np_segs, 2);

    % sq velocity (mm/s = µm/ms)
    calculated(:,4) = sqL./out(:,7);

    % deformed diameter (um): 
    % L_deform, the cell's elongation length in the contraction segment. 
    % Kim et al. 2018 model the deformed cell as an oblate spheroid whose minor axis is set by the contraction width, 
    % V_deform = pi/6*wC*L_deform^2. 
    % volume-equivalent sphere (V_deform = pi/6*d_c^3), 
    % hence:
    % L_deform = sqrt(d_c^3/wC).
    calculated(:,5) = sqrt( ((out(:,4)./out(:,2)*De_c^2*L)./ ...
        (1+0.8*L/De_c*out(:,4)./out(:,2))) ./ wC );

    % wCDI (dimensionless), carry over from original JoVE code
    calculated(:,6) = calculated(:,4)./calculated(:,3) .* calculated(:,1)./ch_height;

    % recovery time (ms)
    calculated(:,7) = out(:,8);

    % recovery rate (per-pulse fit slope, stored in out col 9; using the loop
    % variable fo here took only the last pulse's value and errored outright
    % when the pulse loop ran zero times)
    calculated(:,8) = out(:,9);

    % recovery category
    calculated(:,9) = out(:,12);
    
    % column names & units
    calc_cols = {'diameter', 'strain', 'V_np',          'V_c', 'def_diameter', ...
        'wCDI_original', 'rec_time', 'rec_rate', 'rec_cat'};
    calc_units = {'µm', 'dimensionless', 'mm/s = µm/ms', 'mm/s = µm/ms', 'µm', ...
        'dimensionless', 'ms',      '',         'categorical'};

    %% concatenate pulse data with calculated values
    
    OUT_array = [out(:,1:11), calculated]; % don't include rec_cat twice
    
    % set names & units for the columns
    column_names = [out_cols(1:11), calc_cols];
    column_units = [out_units(1:11), calc_units];
    
end

