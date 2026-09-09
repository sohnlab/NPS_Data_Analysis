function mNPS_showWindow(window_index)
%MNPS_SHOWWINDOW show one candidate window in figure 42.
%   Shades the window's extent on the three signal panels, colors each pulse
%   segment by the role it plays in that window
%       red    = reference (sizing)     yellow = contraction
%       blue   = recovery               grey   = outside this window
%   and draws that window's recovery fit in the bottom panel.
%
%   The handles and per-window fit data are stashed on the figure by
%   mNPS_readData, so only the window index is needed here. Does nothing if
%   figure 42 was never drawn (e.g. an empty pulse).

    figh = findobj('type','figure', 'number',42);
    if isempty(figh) || ~isappdata(figh, 'mNPS_window_data')
        return
    end
    W = getappdata(figh, 'mNPS_window_data');

    C_ref  = [0.85 0.10 0.10]; % red
    C_cont = [0.90 0.70 0.00]; % yellow
    C_rec  = [0.00 0.45 0.74]; % blue
    C_out  = [0.70 0.70 0.70]; % grey

    %% shade the extent of this window on the signal panels

    delete(findobj(figh, 'Tag','mNPS_window_band')); % clear the previous one
    if window_index <= size(W.win_range,1)
        xr = W.win_range(window_index,:);
        for ii = 1:numel(W.ax_sig)
            xregion(W.ax_sig(ii), xr(1), xr(2), 'FaceColor',[0.00 0.60 0.00], ...
                'FaceAlpha',0.10, 'Tag','mNPS_window_band');
        end
    end

    %% color the segments by their role in this window

    roles = [repmat({C_ref},1,W.num_ref), {C_cont}, repmat({C_rec},1,W.num_rec)];
    for ii = 1:size(W.lines,1)
        offset = ii - window_index + 1; % position of this segment in the window
        if offset >= 1 && offset <= numel(roles)
            set(W.lines(ii,:), 'color', roles{offset});
        else
            set(W.lines(ii,:), 'color', C_out);
        end
    end

    %% recovery fit for this window

    if ~isgraphics(W.ax_fit)
        return
    end
    cla(W.ax_fit);
    set(W.ax_fit, 'FontSize',9);
    grid(W.ax_fit, 'on');
    xlabel(W.ax_fit, 'time after contraction (ms)');
    ylabel(W.ax_fit, '\DeltaI');

    if isempty(W.rT) || window_index > size(W.rT,1) || any(isnan(W.rT(window_index,:)))
        title(W.ax_fit, 'recovery fit: not available');
        return
    end

    rT = W.rT(window_index,:);
    rdI = W.rdI(window_index,:);
    p1 = W.p1(window_index);
    p2 = W.p2(window_index);

    hold(W.ax_fit, 'on');
    plot(W.ax_fit, rT, p2+p1*rT, '-', 'color',C_rec, 'linew',1.5);
    plot(W.ax_fit, rT, rdI, 'o', 'color',C_rec, 'markerfacecolor','w', ...
        'markersi',7, 'linew',1.5);
    hold(W.ax_fit, 'off');
    title(W.ax_fit, sprintf('recovery fit, window %d   (slope %+.2e per ms)', ...
        window_index, p1));

    drawnow;

end
