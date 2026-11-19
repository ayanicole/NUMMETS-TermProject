function fig_style(w, h)
% FIG_STYLE  Apply the group's shared figure settings (IEEE column width).
% USAGE : figure; fig_style;           % default 3.5 x 2.6 in
%         figure; fig_style(3.5, 3);   % custom width and height (inches)
    if nargin < 1, w = 3.5; end
    if nargin < 2, h = 2.6; end

    set(gcf, 'Units', 'inches', 'Position', [1 1 w h]);
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 10, 'LineWidth', 0.75, 'Box', 'on', 'TickDir', 'in');
    set(gcf, 'DefaultLineLineWidth', 1.25);
    set(findall(gcf, 'Type', 'line'), 'LineWidth', 1.25);
    grid on;
    colormap(parula);
end
