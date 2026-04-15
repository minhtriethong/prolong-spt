function histogram_dataset_log10(data, xlabel_text, ylabel_text, x_range_min, x_range_max, y_range_min, y_range_max)
% histogram_dataset_log10
% Creates a grid of histogram panels (log10-transformed) with manual axis ranges.
%
% USAGE:
%   histogram_dataset_log10( ...
%       data, xlabel_text, ylabel_text, x_range_min, x_range_max, y_range_min, y_range_max)
%
% INPUTS:
%   data         : Cell array, each cell containing a numeric vector.
%   xlabel_text  : String for X label.
%   ylabel_text  : String for Y label.
%   x_range_min  : Scalar; set NaN/[] to auto from data (after log10).
%   x_range_max  : Scalar; set NaN/[] to auto from data (after log10).
%   y_range_min  : Scalar; set NaN/[] to auto (defaults to 0).
%   y_range_max  : Scalar; set NaN/[] to auto from global histogram peak.

%% Parameters (fixed)
bin_no = 70;
individual_width  = 256;         % width of each histogram panel (pixels)
individual_height = 160;         % height of each histogram panel (pixels)
title_prefix = 'Cell-';          % prefix for each title text
xtick_label_power_of_10 = true;  % display the xtick label in power of 10
apply_log10 = true;              % set to false to disable log10 transformation

%% Input defaults / cleanup
if nargin < 3 || isempty(ylabel_text), ylabel_text = 'Probability density'; end
if nargin < 4 || isempty(x_range_min), x_range_min = NaN; end
if nargin < 5 || isempty(x_range_max), x_range_max = NaN; end
if nargin < 6 || isempty(y_range_min), y_range_min = NaN; end
if nargin < 7 || isempty(y_range_max), y_range_max = NaN; end

% Prepare layout
n = numel(data);

if n > 48
    warning('histogram_dataset_log10:MaxCellsExceeded', ...
        'Input contains %d cells. Only the first 48 cells will be processed.', n);
    data = data(1:48);
    n = 48;
end

if n <= 4
    n_cols = 2;
elseif n <= 9
    n_cols = 3;
elseif n <= 16
    n_cols = 4;
elseif n <= 20
    n_cols = 5;
elseif n <= 30
    n_cols = 6;
elseif n <= 42
    n_cols = 7;
else
    n_cols = 8;
end

n_rows = ceil(n / n_cols);
total_fig_width  = individual_width  * n_cols;
total_fig_height = individual_height * n_rows;
xlabel_text = strrep(xlabel_text, '_', ' ');

%% Optionally apply log10 transformation (and drop non-finite)
if apply_log10
    data = cellfun(@(x) log10(double(x)), data, 'UniformOutput', false);
end
% Remove non-finite values (NaN, Inf, -Inf) from each cell
data = cellfun(@(x) x(isfinite(x)), data, 'UniformOutput', false);

%% Auto-compute x-range across all datasets (after log10), if requested
if isnan(x_range_min) || isnan(x_range_max)
    validMask = ~cellfun(@isempty, data);
    if ~any(validMask)
        error('All elements in <data> are empty after preprocessing.');
    end
    cell_mins = cellfun(@(x) min(x(:)), data(validMask));
    cell_maxs = cellfun(@(x) max(x(:)), data(validMask));

    if isnan(x_range_min), x_range_min = floor(min(cell_mins)); end
    if isnan(x_range_max), x_range_max = ceil(max(cell_maxs));  end
end
% Ensure non-degenerate x-range
if ~isfinite(x_range_min) || ~isfinite(x_range_max)
    error('Computed x-range is not finite. Check input data.');
end
if x_range_max <= x_range_min
    x_range_max = x_range_min + 1;  % minimal span
end

% Use global, fixed bin edges so that y auto-scaling is consistent across panels
bin_edges = linspace(x_range_min, x_range_max, bin_no + 1);

%% Auto-compute y-range (global across panels), if requested
% Use probability-normalized histograms with the same bin edges
if isnan(y_range_min), y_range_min = 0; end
if isnan(y_range_max)
    global_ymax = 0;
    for k = 1:n
        if ~isempty(data{k})
            counts = histcounts(data{k}, 'BinEdges', bin_edges, 'Normalization', 'probability');
            if ~isempty(counts)
                local_max = max(counts);
                if local_max > global_ymax
                    global_ymax = local_max;
                end
            end
        end
    end
    % Round up to 2 decimals for a clean axis; ensure positive span
    y_range_max = max(y_range_min + 0.01, ceil(global_ymax * 100) / 100);
end

%% Figure positioning (avoid top menu overlap)
default_x = 100;
default_y = 100;
offset = 100;
screenSize = get(0, 'ScreenSize');  % [left, bottom, width, height]
screen_height = screenSize(4);
if default_y + total_fig_height > screen_height - offset
    new_y = max(screen_height - offset - total_fig_height, 1);
else
    new_y = default_y;
end

%% Axes margins (normalized) for each panel
panel_outerpos = [0.09, 0.16, 0.98, 0.82];

%% Create overall figure
figure1 = figure('Color',[1 1 1], 'Units','pixels', ...
                 'Position',[default_x, new_y, total_fig_width, total_fig_height]);

%% Loop: Create each histogram in its own uipanel
for i = 1:n
    % Grid position
    col = mod(i-1, n_cols) + 1;
    row = floor((i-1) / n_cols) + 1;

    % Panel pixel position
    panel_left   = (col - 1) * individual_width;
    panel_bottom = total_fig_height - row * individual_height;

    % Panel
    panel = uipanel('Parent', figure1, 'Units', 'pixels', ...
                    'Position', [panel_left, panel_bottom, individual_width, individual_height], ...
                    'BackgroundColor', [1 1 1], 'BorderType', 'none');

    % Axes
    ax = axes('Parent', panel, 'Units', 'normalized', ...
              'OuterPosition', panel_outerpos, ...
              'FontSize', 10, 'TickDir', 'out', 'TickLength', [0.02, 0.025]);
    hold(ax, 'on');

    % Histogram (use fixed bin edges for consistency)
    if ~isempty(data{i})
        histogram(real(double(data{i})), 'Parent', ax, ...
                  'FaceColor', [1 1 1], 'FaceAlpha', 0, ...
                  'Normalization', 'probability', ...
                  'BinEdges', bin_edges);
    else
        % Plot an empty histogram to keep layout consistent
        histogram([], 'Parent', ax, 'Normalization', 'probability', 'BinEdges', bin_edges, ...
                  'FaceColor', [1 1 1], 'FaceAlpha', 0);
    end

    % Axis limits
    xlim(ax, [x_range_min, x_range_max]);
    ylim(ax, [y_range_min, y_range_max]);

    % Optional power-of-10 x tick labels (remember: x-axis is log10 scale values)
    if xtick_label_power_of_10
        xt = get(ax, 'XTick');
        labels = arrayfun(@(x) sprintf('10^{%g}', x), xt, 'UniformOutput', false);
        set(ax, 'XTickLabel', labels);
    end

    % Labels
    xlabel_handle = xlabel(ax, xlabel_text, 'FontSize', 10);
    set(xlabel_handle, 'Units', 'normalized', 'Position', [0.5, -0.35, 0]);

    ylabel_handle = ylabel(ax, ylabel_text, 'FontSize', 10);
    set(ylabel_handle, 'Units', 'normalized', 'Position', [-0.18, 0.5, 0]);

    % Title
    title_text = sprintf('%s%02d', title_prefix, i);
    title_handle = title(ax, title_text, 'FontSize', 10, ...
                         'FontWeight', 'normal', 'HorizontalAlignment', 'left');
    set(title_handle, 'Units', 'normalized', 'Position', [0.04, 0.85, 0]);

    hold(ax, 'off');
end

%% Optional: fill unused grid cells with blank panels
total_slots = n_rows * n_cols;
if n < total_slots
    for i = n+1 : total_slots
        col = mod(i-1, n_cols) + 1;
        row = floor((i-1) / n_cols) + 1;
        panel_left   = (col - 1) * individual_width;
        panel_bottom = total_fig_height - row * individual_height;
        uipanel('Parent', figure1, 'Units', 'pixels', ...
                'Position', [panel_left, panel_bottom, individual_width, individual_height], ...
                'BackgroundColor', [1 1 1], 'BorderType', 'none');
    end
end

end