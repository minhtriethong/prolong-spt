function trackingMovie(processingInfo)
% Fast renderer: builds RGB frames in memory and writes directly to TIFF.
% Output: single movie with raw frame + yellow trajectories + yellow outline circles.
%
% Frames are read directly from processingInfo.sourceImagePath using Bio-Formats.

%% ---------- Config ----------
trail_len       = 50;      % window length in frames
line_width      = 1;       % pixels for lines and circle outlines
dot_radius      = 3;       % circle radius in pixels
write_every     = 1;       % 1 = every frame
use_big_tiff    = false;   % set true if >4GB output
gap_fill_max    = 5;       % default gap fill max
gap_fill_method = 'linear';
ring_alpha      = 0.5;     % 0 = fully transparent, 1 = fully opaque
trail_end_alpha = 0.0;     % oldest trail segment transparency
%USE_CVT         = (exist('insertShape','file') == 2);
USE_CVT         = false;
YEL = uint8([255 255 0]);  % single yellow color for everything

%% ---------- processingInfo ----------
tracksCoordinates = processingInfo.tracksCoordinates;
sourceImagePath   = char(processingInfo.sourceImagePath);
imageName         = char(processingInfo.imageName);
T_movie           = double(processingInfo.sizeT);
H                 = double(processingInfo.sizeY);
W                 = double(processingInfo.sizeX);
cmax              = double(processingInfo.AppInputSettings.ContrastMaxVal);
cmin              = double(processingInfo.AppInputSettings.ContrastMinVal);

% Optional override from processingInfo.AppInputSettings.RingAlpha
% if isfield(processingInfo, 'AppInputSettings') && ...
%         isfield(processingInfo.AppInputSettings, 'RingAlpha') && ...
%         ~isempty(processingInfo.AppInputSettings.RingAlpha)
%     ring_alpha = double(processingInfo.AppInputSettings.RingAlpha);
% end
ring_alpha = max(0, min(1, ring_alpha));
trail_end_alpha = max(0, min(1, trail_end_alpha));

%% ---------- Paths / output ----------
[source_folder, ~, ~] = fileparts(sourceImagePath);
out_folder = fullfile(source_folder, 'spt-check');
if ~exist(out_folder, 'dir')
    mkdir(out_folder);
end

out_name = sprintf('%s_SPT.tif', imageName);
tiff_filename = fullfile(out_folder, out_name);
if exist(tiff_filename, 'file')
    delete(tiff_filename);
end

%% ---------- Bio-Formats reader ----------
if exist('bfGetReader', 'file') ~= 2 || exist('bfGetPlane', 'file') ~= 2
    error('Bio-Formats MATLAB helpers bfGetReader and bfGetPlane must be on the MATLAB path.');
end

reader = bfGetReader(sourceImagePath);
readerCleanup = onCleanup(@() reader.close()); %#ok<NASGU>

seriesIdx = 1;
if isfield(processingInfo, 'seriesIndex') && ~isempty(processingInfo.seriesIndex)
    seriesIdx = double(processingInfo.seriesIndex);
end
reader.setSeries(seriesIdx - 1);

zIdx = 1;
if isfield(processingInfo, 'zIndex') && ~isempty(processingInfo.zIndex)
    zIdx = double(processingInfo.zIndex);
end

cIdx = 1;
if isfield(processingInfo, 'channelIndex') && ~isempty(processingInfo.channelIndex)
    cIdx = double(processingInfo.channelIndex);
end

readerT = double(reader.getSizeT());
readerZ = double(reader.getSizeZ());
readerC = double(reader.getSizeC());

if zIdx < 1 || zIdx > readerZ
    error('Requested zIndex (%d) is outside sourceImagePath Z range [1, %d].', zIdx, readerZ);
end
if cIdx < 1 || cIdx > readerC
    error('Requested channelIndex (%d) is outside sourceImagePath C range [1, %d].', cIdx, readerC);
end
if T_movie < 1 || T_movie > readerT
    error('processingInfo.sizeT (%d) is outside sourceImagePath T range [1, %d].', T_movie, readerT);
end

%% ---------- Gap filling ----------
[tracksFilled, ~] = single_particle_tracking_gap_filling_v01( ...
    tracksCoordinates, gap_fill_max, gap_fill_method);

%% ---------- Data ----------
X = tracksFilled.X;   % [T x M]
Y = tracksFilled.Y;   % [T x M]

if isempty(X)
    error('processingInfo.tracksCoordinates.X is empty.');
end
if ~isequal(size(X), size(Y))
    error('processingInfo.tracksCoordinates.X and .Y must have the same size.');
end

[Txy, ~] = size(X);

%% ---------- Render loop ----------
firstWritten = false;

for i = 1:write_every:T_movie
    % ----- Load + contrast-adjust the raw frame -----
    img_plane = read_bioformats_plane(reader, i, zIdx, cIdx);

    if size(img_plane, 1) ~= H || size(img_plane, 2) ~= W
        error('Frame %d size [%d x %d] does not match processingInfo size [%d x %d].', ...
            i, size(img_plane,1), size(img_plane,2), H, W);
    end

    img_u8 = apply_processinginfo_contrast(img_plane, cmin, cmax);

    % Single output frame base (RGB)
    frame_rgb = repmat(img_u8, [1 1 3]);

    % Clamp to tracking axis
    ti     = min(i, Txy);
    wStart = max(1, ti - trail_len + 1);
    win    = wStart:ti;
    Lw     = numel(win);
    sc     = max(0, Lw - 1);  % #segments in window

    % ----- Draw trajectory segments (yellow, fading) -----
    % Draw oldest -> newest so newest stays on top.
    % Oldest segment alpha = trail_end_alpha
    % Newest segment alpha = ring_alpha
    if sc > 0
        if sc == 1
            alphaRank = ring_alpha;
        else
            alphaRank = linspace(trail_end_alpha, ring_alpha, sc);  % oldest -> newest
        end

        for k = 1:sc
            x1 = X(win(k),   :);
            y1 = Y(win(k),   :);
            x2 = X(win(k+1), :);
            y2 = Y(win(k+1), :);

            valid = ~(isnan(x1) | isnan(y1) | isnan(x2) | isnan(y2));
            if any(valid)
                lines = [x1(valid)' y1(valid)' x2(valid)' y2(valid)'];  % Nx4
                frame_rgb = draw_lines_masked(frame_rgb, lines, double(YEL), alphaRank(k), line_width);
            end
        end
    end

    % ----- Draw current dots as outline circles (yellow) -----
    xcur = X(ti, :);
    ycur = Y(ti, :);
    mask = ~(isnan(xcur) | isnan(ycur));
    if any(mask)
        centers = [xcur(mask)' ycur(mask)' repmat(dot_radius, [nnz(mask), 1])];  % Nx3 [x y r]
        frame_rgb = draw_circles_outline_masked(frame_rgb, centers, double(YEL), ring_alpha, line_width);
    end

    % ----- Write to TIFF -----
    if ~firstWritten
        if use_big_tiff
            imwrite(frame_rgb, tiff_filename, 'tif', ...
                'Compression', 'none', 'WriteMode', 'overwrite', 'BigTIFF', true);
        else
            imwrite(frame_rgb, tiff_filename, 'tif', ...
                'Compression', 'none', 'WriteMode', 'overwrite');
        end
        firstWritten = true;
    else
        imwrite(frame_rgb, tiff_filename, 'tif', ...
            'Compression', 'none', 'WriteMode', 'append');
    end
end
end

%% -------------------- Bio-Formats helpers --------------------

function img = read_bioformats_plane(reader, tIdx, zIdx, cIdx)
% Bio-Formats getIndex uses 0-based Z/C/T.
planeIdx = reader.getIndex(zIdx - 1, cIdx - 1, tIdx - 1) + 1;
img = bfGetPlane(reader, planeIdx);

if ndims(img) == 3
    if size(img, 3) == 1
        img = img(:, :, 1);
    else
        img = rgb2gray(img);
    end
end
end

function img_u8 = apply_processinginfo_contrast(img, cmin, cmax)
% Uses processingInfo.AppInputSettings.ContrastMinVal / ContrastMaxVal.
% Supports either normalized [0,1] limits or raw-value limits.

if ndims(img) == 3
    if size(img, 3) == 1
        img = img(:, :, 1);
    else
        img = rgb2gray(img);
    end
end

cmin = double(cmin);
cmax = double(cmax);

if cmax <= cmin
    error('processingInfo.AppInputSettings.ContrastMaxVal must be greater than ContrastMinVal.');
end

if isinteger(img)
    img01 = im2double(img);

    if cmin >= 0 && cmax <= 1
        low  = cmin;
        high = cmax;
    else
        classMax = double(intmax(class(img)));
        low  = cmin / classMax;
        high = cmax / classMax;
    end

    low  = max(0, min(1, low));
    high = max(0, min(1, high));

    if high <= low
        error('Invalid contrast window after scaling ContrastMinVal / ContrastMaxVal.');
    end

    img_u8 = im2uint8(imadjust(img01, [low high], [0 1]));
else
    img = double(img);

    if min(img(:)) >= 0 && max(img(:)) <= 1 && cmin >= 0 && cmax <= 1
        low  = max(0, min(1, cmin));
        high = max(0, min(1, cmax));

        if high <= low
            error('Invalid normalized contrast window.');
        end

        img_u8 = im2uint8(imadjust(img, [low high], [0 1]));
    else
        img01 = (img - cmin) / (cmax - cmin);
        img01 = min(max(img01, 0), 1);
        img_u8 = im2uint8(img01);
    end
end
end

%% -------------------- Drawing helpers --------------------

function rgb = draw_lines_masked(rgb, lines, color, alpha, line_width)
% Rasterize all line pixels into a mask, then alpha-blend once per pixel.
[H, W, ~] = size(rgb);
if isempty(lines) || alpha <= 0
    return;
end

mask = false(H, W);
brush = max(1, round(line_width));
halfB = floor((brush - 1) / 2);

for i = 1:size(lines, 1)
    x1 = lines(i, 1);
    y1 = lines(i, 2);
    x2 = lines(i, 3);
    y2 = lines(i, 4);

    nS = max(1, round(max(abs(x2 - x1), abs(y2 - y1)))) + 1;
    xs = round(linspace(x1, x2, nS));
    ys = round(linspace(y1, y2, nS));

    in = xs >= 1 & xs <= W & ys >= 1 & ys <= H;
    xs = xs(in);
    ys = ys(in);

    if isempty(xs)
        continue;
    end

    pts = unique([ys(:) xs(:)], 'rows', 'stable');
    ys = pts(:,1);
    xs = pts(:,2);

    for k = 1:numel(xs)
        r1 = max(1, ys(k) - halfB);
        r2 = min(H, ys(k) + halfB);
        c1 = max(1, xs(k) - halfB);
        c2 = min(W, xs(k) + halfB);
        mask(r1:r2, c1:c2) = true;
    end
end

rgb = blend_mask(rgb, mask, color, alpha);
end

function rgb = draw_circles_outline_masked(rgb, centers, color, alpha, line_width)
% Rasterize circle outlines into a mask, then alpha-blend once per pixel.
[H, W, ~] = size(rgb);
if isempty(centers) || alpha <= 0
    return;
end

mask = false(H, W);
brush = max(1, round(line_width));
halfB = floor((brush - 1) / 2);

for i = 1:size(centers, 1)
    cx  = centers(i, 1);
    cy  = centers(i, 2);
    rad = centers(i, 3);

    % Oversample circumference to avoid gaps
    nPts = max(32, round(4 * pi * rad));
    th   = linspace(0, 2*pi, nPts + 1);

    xs = round(cx + rad * cos(th));
    ys = round(cy + rad * sin(th));

    in = xs >= 1 & xs <= W & ys >= 1 & ys <= H;
    xs = xs(in);
    ys = ys(in);

    if isempty(xs)
        continue;
    end

    pts = unique([ys(:) xs(:)], 'rows', 'stable');
    ys = pts(:,1);
    xs = pts(:,2);

    for k = 1:numel(xs)
        r1 = max(1, ys(k) - halfB);
        r2 = min(H, ys(k) + halfB);
        c1 = max(1, xs(k) - halfB);
        c2 = min(W, xs(k) + halfB);
        mask(r1:r2, c1:c2) = true;
    end
end

rgb = blend_mask(rgb, mask, color, alpha);
end

function rgb = blend_mask(rgb, mask, color, alpha)
% Alpha-blend one color into all pixels in mask exactly once.
if ~any(mask(:)) || alpha <= 0
    return;
end

alpha = max(0, min(1, double(alpha)));
color = double(color(:))';  % 1x3

for ch = 1:3
    plane = double(rgb(:,:,ch));
    plane(mask) = alpha * color(ch) + (1 - alpha) * plane(mask);
    rgb(:,:,ch) = uint8(plane);
end
end

%% -------------------- Gap-filling helpers --------------------

function [tracksFilled, filledMask] = single_particle_tracking_gap_filling_v01(tracksCoordinates, maxGap, method)
% Fills short (<= maxGap) interior NaN runs in tracksCoordinates.X and .Y.
% - Interpolates along time (rows) independently per trajectory (column).
% - Leaves leading/trailing NaNs and long gaps (>maxGap) as NaN.
% - Returns the filled tracks and logical masks of where values were imputed.
%
% Inputs:
%   tracksCoordinates : struct with fields X [T x M], Y [T x M] (I is passed through)
%   maxGap            : (optional) max consecutive NaNs to fill (default 5)
%   method            : (optional) 'linear' (default) | 'pchip' | 'spline'
%
% Outputs:
%   tracksFilled : struct like input, with X/Y gaps filled where eligible
%   filledMask   : struct with fields X/Y (logical [T x M]) of filled samples

if nargin < 2 || isempty(maxGap), maxGap = 5; end
if nargin < 3 || isempty(method), method = 'linear'; end

X = tracksCoordinates.X;
Y = tracksCoordinates.Y;

% Prefer fillmissing if available (fast, robust); otherwise manual fallback
if exist('fillmissing', 'file') == 2
    Xf = fillmissing(X, method, 1, 'MaxGap', maxGap);
    Yf = fillmissing(Y, method, 1, 'MaxGap', maxGap);
else
    Xf = fill_short_gaps_manual(X, maxGap, method);
    Yf = fill_short_gaps_manual(Y, maxGap, method);
end

filledMask.X = isnan(X) & ~isnan(Xf);
filledMask.Y = isnan(Y) & ~isnan(Yf);

tracksFilled   = tracksCoordinates;
tracksFilled.X = Xf;
tracksFilled.Y = Yf;

if isfield(tracksCoordinates, 'I')
    tracksFilled.I = tracksCoordinates.I;
end
end

function Zf = fill_short_gaps_manual(Z, maxGap, method)
% Fill <= maxGap NaN runs inside the series (no edges), per column.
[T, M] = size(Z);
Zf = Z;

for m = 1:M
    z = Z(:, m);
    isn = isnan(z);

    if ~any(isn)
        Zf(:, m) = z;
        continue;
    end

    % Find NaN runs
    d = diff([0; isn; 0]);
    starts = find(d == 1);
    ends   = find(d == -1) - 1;
    lens   = ends - starts + 1;

    for k = 1:numel(starts)
        s = starts(k);
        e = ends(k);
        L = lens(k);

        % Only fill if strictly interior and short
        if L <= maxGap && s > 1 && e < T && ~isnan(z(s-1)) && ~isnan(z(e+1))
            t_known = [s-1; e+1];
            y_known = [z(s-1); z(e+1)];
            t_query = (s:e).';

            z(t_query) = interp1(t_known, y_known, t_query, method);
        end
    end

    Zf(:, m) = z;
end
end