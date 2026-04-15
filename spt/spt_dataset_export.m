function spt_dataset_export(Results_folderpath, immobile_threshold)
%spt_dataset_export Build, export, and plot an SPT dataset summary from a Results folder.
%
%   Usage
%   -----
%   spt_dataset_export(Results_folderpath)
%   spt_dataset_export(Results_folderpath, immobile_threshold)
%
%   Input
%   -----
%   Results_folderpath : character vector or string scalar
%       Path to a Results folder. Each immediate subfolder inside this
%       folder is treated as one dataset entry and should contain a file
%       named processingInfo.mat.
%
%   immobile_threshold : numeric scalar, optional
%       Threshold applied as:
%           log10(processingInfo.track_features.D.D_track) < immobile_threshold
%       If provided, dataset.mean.immobilefraction is calculated and
%       included in datasetInfo.mat and datasetInfo.csv.
%
%   Extra outputs written by this version
%   -------------------------------------
%   In each subfolder:
%       traj_map.fig
%
%   In Results_folderpath:
%       traj_map.png   % combined flattened grid of trajectory maps
%
%   Required external functions
%   ---------------------------
%   histogram_dataset_log10
%   traj_map_plot

if nargin < 1 || nargin > 2
    error('spt_dataset_export:InvalidInputCount', ...
        'Usage: spt_dataset_export(Results_folderpath) or spt_dataset_export(Results_folderpath, immobile_threshold).');
end

if ~(ischar(Results_folderpath) || isstring(Results_folderpath))
    error('spt_dataset_export:InvalidInput', ...
        'Results_folderpath must be a character vector or string scalar.');
end

Results_folderpath = char(Results_folderpath);

if ~isfolder(Results_folderpath)
    error('spt_dataset_export:FolderNotFound', ...
        'Results folder does not exist: %s', Results_folderpath);
end

calculateImmobileFraction = (nargin >= 2) && ~isempty(immobile_threshold);

if calculateImmobileFraction
    if ~(isnumeric(immobile_threshold) && isscalar(immobile_threshold) && isfinite(immobile_threshold))
        error('spt_dataset_export:InvalidImmobileThreshold', ...
            'immobile_threshold must be a finite numeric scalar when provided.');
    end
end

mu = char(956);
sup2 = char(178);

header_D = ['D (' mu 'm' sup2 '/s)'];
header_travelDistances = ['travelDistances (' mu 'm)'];
header_travelTime = 'travelTime (s)';
header_majorAxis = ['majorAxis (' mu 'm)'];
header_minorAxis = ['minorAxis (' mu 'm)'];
header_trackNo = 'trackNo';
header_immobileFraction = 'immobile fraction';

label_D = ['Diffusion coefficient (' mu 'm' sup2 '/s)'];
label_majorAxis = ['Major axis (' mu 'm)'];
label_minorAxis = ['Minor axis (' mu 'm)'];
label_travelDistance = ['Travel distance (' mu 'm)'];
label_travelTime = 'Travel time (s)';
label_y = 'Probability density';

% Plot settings (edit if needed)
plotSettings.D.x = [-7, 4];
plotSettings.D.y = [0, 0.2];

plotSettings.majorAxis.x = [-3, 1];
plotSettings.majorAxis.y = [0, 0.15];

plotSettings.minorAxis.x = [-4, 1];
plotSettings.minorAxis.y = [0, 0.15];

plotSettings.travelDistances.x = [-3, 3];
plotSettings.travelDistances.y = [0, 0.15];

plotSettings.travelTime.x = [-1, 2];
plotSettings.travelTime.y = [0, 0.4];

% Combined trajectory map settings
maxTrajGridSubfolders = 50;
trajGridPaddingPx = 12;
trajGridMaxBytes = 50 * 1024 * 1024;   % 50 MB target limit when trajGridQuality <= 1
trajGridQuality = 1.0;                 % 1 = default.
                                       % >1 = higher quality traj_map.png and ignore 50 MB cap.
                                       % <1 = smaller/lighter combined traj_map.png.
trajTileAspect = 0.78;                 % height / width
trajGridRenderDPI = 500;               % internal render DPI for flattened PNG tiles

validateattributes(trajGridQuality, {'numeric'}, ...
    {'scalar','real','finite','>',0}, mfilename, 'trajGridQuality');

trajTileMinWidthPx = max(250, round(500 * trajGridQuality));
trajTileMaxWidthPx = max(trajTileMinWidthPx, round(900 * trajGridQuality));
trajGridMaxWidthPx = max(2000, round(7000 * trajGridQuality));

if trajGridQuality > 1
    trajGridEffectiveMaxBytes = inf;   % user explicitly asked for higher quality over size cap
else
    trajGridEffectiveMaxBytes = trajGridMaxBytes;
end

% Get immediate subfolders only
dirInfo = dir(Results_folderpath);
isSubfolder = [dirInfo.isdir] & ~ismember({dirInfo.name}, {'.', '..'});
subfolders = dirInfo(isSubfolder);

if isempty(subfolders)
    error('spt_dataset_export:NoSubfolders', ...
        'No subfolders found in: %s', Results_folderpath);
end

nFolders = numel(subfolders);

if nFolders > maxTrajGridSubfolders
    warning('spt_dataset_export:TooManySubfoldersForTrajGrid', ...
        ['Found %d subfolders. The combined traj_map.png will include only ' ...
         'the first %d subfolders in folder order.'], ...
        nFolders, maxTrajGridSubfolders);
end

% Precompute tile size for the combined grid
previewCount = min(nFolders, maxTrajGridSubfolders);
previewNcol = pick_traj_grid_ncol(previewCount);
trajTileWidthPx = floor((trajGridMaxWidthPx - trajGridPaddingPx * (previewNcol + 1)) / previewNcol);
trajTileWidthPx = min(trajTileMaxWidthPx, max(trajTileMinWidthPx, trajTileWidthPx));
trajTileHeightPx = max(round(380 * trajGridQuality), round(trajTileWidthPx * trajTileAspect));
trajFigSizePx = [trajTileWidthPx, trajTileHeightPx];

dataset = struct();

dataset.raw = struct();
dataset.raw.D = cell(nFolders, 1);
dataset.raw.travelDistances = cell(nFolders, 1);
dataset.raw.travelTime = cell(nFolders, 1);
dataset.raw.majorAxis = cell(nFolders, 1);
dataset.raw.minorAxis = cell(nFolders, 1);

dataset.mean = struct();
dataset.mean.D = nan(nFolders, 1);
dataset.mean.travelDistances = nan(nFolders, 1);
dataset.mean.travelTime = nan(nFolders, 1);
dataset.mean.majorAxis = nan(nFolders, 1);
dataset.mean.minorAxis = nan(nFolders, 1);
dataset.mean.trackNo = nan(nFolders, 1);

if calculateImmobileFraction
    dataset.mean.immobilefraction = nan(nFolders, 1);
end

hasTrajMapFunction = exist('traj_map_plot', 'file') == 2;
if ~hasTrajMapFunction
    warning('spt_dataset_export:TrajMapFunctionMissing', ...
        ['traj_map_plot was not found on the MATLAB path. ' ...
         'Subfolder traj_map.fig export and combined traj_map.png will be skipped.']);
end

trajMapFlatImages = cell(0,1);

for i = 1:nFolders
    folderName = subfolders(i).name;
    subfolderPath = fullfile(Results_folderpath, folderName);
    matFile = fullfile(subfolderPath, 'processingInfo.mat');

    if ~isfile(matFile)
        warning('spt_dataset_export:MissingFile', ...
            'File not found, skipped: %s', matFile);
        continue;
    end

    try
        S = load(matFile, 'processingInfo');

        if ~isfield(S, 'processingInfo')
            warning('spt_dataset_export:MissingVariable', ...
                'Variable "processingInfo" not found in: %s', matFile);
            continue;
        end

        processingInfo = S.processingInfo;

        % Save visible traj_map.fig in each subfolder, then flatten for contact sheet
        if hasTrajMapFunction
            try
                flatImg = generate_traj_map_fig( ...
                    processingInfo, subfolderPath, 'traj_map', ...
                    trajFigSizePx, trajGridRenderDPI);

                if i <= maxTrajGridSubfolders && ~isempty(flatImg)
                    trajMapFlatImages{end+1,1} = flatImg; %#ok<AGROW>
                end
            catch MEtraj
                warning('spt_dataset_export:TrajMapFailed', ...
                    'Failed to create trajectory map for %s\n%s', matFile, MEtraj.message);
            end
        end

        if ~isfield(processingInfo, 'track_features') || isempty(processingInfo.track_features)
            warning('spt_dataset_export:MissingTrackFeatures', ...
                'Field "track_features" not found or empty in: %s', matFile);
            continue;
        end

        track_features = processingInfo.track_features;

        dataset.raw.D{i,1} = track_features.D.D_track_rolling_combined;
        dataset.raw.travelDistances{i,1} = track_features.travelDistances;
        dataset.raw.travelTime{i,1} = track_features.travelTime;
        dataset.raw.majorAxis{i,1} = track_features.majorAxis;
        dataset.raw.minorAxis{i,1} = track_features.minorAxis;

        dataset.mean.D(i,1) = double(track_features.D.D_cell);
        dataset.mean.travelDistances(i,1) = mean(track_features.travelDistances(:));
        dataset.mean.travelTime(i,1) = mean(track_features.travelTime(:));
        dataset.mean.majorAxis(i,1) = mean(track_features.majorAxis(:));
        dataset.mean.minorAxis(i,1) = mean(track_features.minorAxis(:));
        dataset.mean.trackNo(i,1) = length(processingInfo.track_features.D.D_track);

        if calculateImmobileFraction
            D_track = processingInfo.track_features.D.D_track;

            if isempty(D_track)
                dataset.mean.immobilefraction(i,1) = NaN;
            else
                dataset.mean.immobilefraction(i,1) = ...
                    length(find(log10(D_track) < immobile_threshold)) / length(D_track); %#ok<FNDSB>
            end
        end

    catch ME
        warning('spt_dataset_export:ProcessingFailed', ...
            'Failed to process %s\n%s', matFile, ME.message);
    end
end

save(fullfile(Results_folderpath, 'datasetInfo.mat'), 'dataset');

headers = { ...
    header_D, ...
    header_travelDistances, ...
    header_travelTime, ...
    header_majorAxis, ...
    header_minorAxis, ...
    header_trackNo};

summaryMatrix = [ ...
    dataset.mean.D(:,1), ...
    dataset.mean.travelDistances(:,1), ...
    dataset.mean.travelTime(:,1), ...
    dataset.mean.majorAxis(:,1), ...
    dataset.mean.minorAxis(:,1), ...
    dataset.mean.trackNo(:,1)];

if calculateImmobileFraction
    headers{end+1} = header_immobileFraction;
    summaryMatrix = [summaryMatrix, dataset.mean.immobilefraction(:,1)];
end

csvFile = fullfile(Results_folderpath, 'datasetInfo.csv');
fid = fopen(csvFile, 'w', 'n', 'UTF-8');

if fid == -1
    error('spt_dataset_export:CSVOpenFailed', ...
        'Cannot open CSV file for writing: %s', csvFile);
end

cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>

fwrite(fid, uint8([239 187 191]), 'uint8');

headerFormat = [repmat('%s,', 1, numel(headers)-1), '%s\r\n'];
fprintf(fid, headerFormat, headers{:});

rowFormat = [repmat('%.15g,', 1, size(summaryMatrix, 2)-1), '%.15g\r\n'];
for i = 1:size(summaryMatrix, 1)
    fprintf(fid, rowFormat, summaryMatrix(i, :));
end

% Build combined trajectory-map grid
if hasTrajMapFunction && ~isempty(trajMapFlatImages)
    try
        generate_traj_map_grid( ...
            trajMapFlatImages, Results_folderpath, 'traj_map', ...
            trajGridPaddingPx, trajGridEffectiveMaxBytes);
    catch MEgrid
        warning('spt_dataset_export:TrajGridFailed', ...
            'Combined traj_map.png could not be created.\n%s', MEgrid.message);
    end
end

hasHistogramFunction = exist('histogram_dataset_log10', 'file') ~= 0;
if ~hasHistogramFunction
    warning('spt_dataset_export:PlotFunctionMissing', ...
        ['datasetInfo.mat and datasetInfo.csv were saved, but histogram plot files ' ...
         'were skipped because histogram_dataset_log10 was not found on the MATLAB path.']);
    return;
end

generate_dataset_plot(dataset.raw.D, label_D, label_y, ...
    plotSettings.D.x, plotSettings.D.y, Results_folderpath, 'D_hist');

generate_dataset_plot(dataset.raw.majorAxis, label_majorAxis, label_y, ...
    plotSettings.majorAxis.x, plotSettings.majorAxis.y, Results_folderpath, 'Major_axis');

generate_dataset_plot(dataset.raw.minorAxis, label_minorAxis, label_y, ...
    plotSettings.minorAxis.x, plotSettings.minorAxis.y, Results_folderpath, 'Minor_axis');

generate_dataset_plot(dataset.raw.travelDistances, label_travelDistance, label_y, ...
    plotSettings.travelDistances.x, plotSettings.travelDistances.y, Results_folderpath, 'Travel_distance');

generate_dataset_plot(dataset.raw.travelTime, label_travelTime, label_y, ...
    plotSettings.travelTime.x, plotSettings.travelTime.y, Results_folderpath, 'Travel_time');

end

function generate_dataset_plot(dataCell, xLabelText, yLabelText, xRange, yRange, outputFolder, baseName)
histogram_dataset_log10( ...
    dataCell, xLabelText, yLabelText, ...
    xRange(1), xRange(2), yRange(1), yRange(2));

fig = gcf;
drawnow;

save_plot(fullfile(outputFolder, [baseName '.jpg']));
saveas(gcf, fullfile(outputFolder, [baseName '.pdf']));

close(fig);
end

function flatImg = generate_traj_map_fig(processingInfo, outputFolder, baseName, figSizePx, renderDPI)
% Create visible trajectory-map FIG and return a flattened RGB image.

figFile = fullfile(outputFolder, [baseName '.fig']);
flatImg = [];

hFig = traj_map_plot(processingInfo);
figCleaner = onCleanup(@() close_if_valid(hFig));

% Make figure visible and standardize size before saving/rendering
% try
%     set(hFig, 'Visible', 'on', 'Units', 'pixels');
%     pos = get(hFig, 'Position');
%     pos(3) = figSizePx(1);
%     pos(4) = figSizePx(2);
%     set(hFig, 'Position', pos);
% catch
% end

drawnow;

try
    savefig(hFig, figFile, 'compact');
catch
    savefig(hFig, figFile);
end

flatImg = render_figure_to_rgb(hFig, figSizePx, renderDPI);
end

function img = render_figure_to_rgb(hFig, renderSizePx, renderDPI)
% Render a figure to an RGB image at a controlled output size.

tmpPng = [tempname '.png'];
tmpCleaner = onCleanup(@() delete_if_exists(tmpPng)); %#ok<NASGU>

img = [];

try
    % set(hFig, 'Units', 'pixels');
    % pos = get(hFig, 'Position');
    % pos(3) = renderSizePx(1);
    % pos(4) = renderSizePx(2);
    % set(hFig, 'Position', pos);

    widthIn  = renderSizePx(1) / renderDPI;
    heightIn = renderSizePx(2) / renderDPI;

    set(hFig, 'PaperUnits', 'inches');
    set(hFig, 'PaperPositionMode', 'manual');
    set(hFig, 'PaperPosition', [0 0 widthIn heightIn]);
    set(hFig, 'PaperSize', [widthIn heightIn]);

    drawnow;

    print(hFig, tmpPng, '-dpng', sprintf('-r%d', renderDPI));
    img = imread(tmpPng);

catch
    % Fallback
    drawnow;
    fr = getframe(hFig);
    img = frame2im(fr);
end

img = force_rgb_uint8(img);
end

function generate_traj_map_grid(flatImages, outputFolder, baseName, paddingPx, maxBytes)
% Combine flattened trajectory-map images into a left-to-right, top-to-bottom grid.

n = numel(flatImages);
if n < 1
    return;
end

ncol = pick_traj_grid_ncol(n);
nrow = ceil(n / ncol);

flatImages = cellfun(@force_rgb_uint8, flatImages, 'UniformOutput', false);

tileHeights = cellfun(@(I) size(I,1), flatImages);
tileWidths  = cellfun(@(I) size(I,2), flatImages);

tileH = max(tileHeights);
tileW = max(tileWidths);

canvasH = paddingPx + nrow * tileH + (nrow - 1) * paddingPx + paddingPx;
canvasW = paddingPx + ncol * tileW + (ncol - 1) * paddingPx + paddingPx;

canvas = uint8(255 * ones(canvasH, canvasW, 3, 'uint8'));

for k = 1:n
    row = floor((k - 1) / ncol) + 1;
    col = mod(k - 1, ncol) + 1;

    tile = flatImages{k};
    th = size(tile, 1);
    tw = size(tile, 2);

    r0 = paddingPx + (row - 1) * (tileH + paddingPx) + 1;
    c0 = paddingPx + (col - 1) * (tileW + paddingPx) + 1;

    rr = r0 + floor((tileH - th) / 2);
    cc = c0 + floor((tileW - tw) / 2);

    canvas(rr:rr+th-1, cc:cc+tw-1, :) = tile;
end

outFile = fullfile(outputFolder, [baseName '.png']);
save_png_under_size_limit(canvas, outFile, maxBytes);
end

function ncol = pick_traj_grid_ncol(n)
% User-specified mapping:
% <=2  -> 2
% <=6  -> 3
% <=12 -> 4
% <=15 -> 5
% <=24 -> 6
% <=28 -> 7
% <=40 -> 8
% <=45 -> 9
% <=50 -> 10

if n <= 2
    ncol = 2;
elseif n <= 6
    ncol = 3;
elseif n <= 12
    ncol = 4;
elseif n <= 15
    ncol = 5;
elseif n <= 24
    ncol = 6;
elseif n <= 28
    ncol = 7;
elseif n <= 40
    ncol = 8;
elseif n <= 45
    ncol = 9;
else
    ncol = 10;
end
end

function save_png_under_size_limit(img, outFile, maxBytes)
% Save a PNG and downsample if needed to keep the file size under maxBytes.
% If maxBytes is Inf, save at full size without downsampling.

imgWork = force_rgb_uint8(img);

imwrite(imgWork, outFile, 'png');
d = dir(outFile);

if isempty(d)
    error('spt_dataset_export:PNGSaveFailed', ...
        'Combined trajectory-map PNG could not be written: %s', outFile);
end

if isinf(maxBytes)
    return;
end

wasDownsampled = false;

while d.bytes > maxBytes && min(size(imgWork,1), size(imgWork,2)) > 400
    wasDownsampled = true;

    scale = sqrt(double(maxBytes) / double(d.bytes)) * 0.95;
    scale = max(min(scale, 0.95), 0.60);

    newH = max(1, floor(size(imgWork,1) * scale));
    newW = max(1, floor(size(imgWork,2) * scale));

    imgWork = downsample_rgb_image(imgWork, newH, newW);

    imwrite(imgWork, outFile, 'png');
    d = dir(outFile);

    if isempty(d)
        error('spt_dataset_export:PNGSaveFailed', ...
            'Combined trajectory-map PNG could not be written: %s', outFile);
    end
end

if wasDownsampled
    warning('spt_dataset_export:TrajGridDownsampled', ...
        ['Combined traj_map.png was downsampled automatically to keep the ' ...
         'file size under %.1f MB.'], maxBytes / 1024 / 1024);
end
end

function imgOut = downsample_rgb_image(imgIn, newH, newW)
% Downsample RGB image. Uses imresize when available, otherwise indexing fallback.

if exist('imresize', 'file') == 2
    try
        imgOut = imresize(imgIn, [newH newW], 'bicubic');
        imgOut = force_rgb_uint8(imgOut);
        return;
    catch
    end
end

rowIdx = round(linspace(1, size(imgIn,1), newH));
colIdx = round(linspace(1, size(imgIn,2), newW));
imgOut = imgIn(rowIdx, colIdx, :);
end

function img = force_rgb_uint8(img)
% Ensure output is uint8 RGB.

if isempty(img)
    img = uint8([]);
    return;
end

if isa(img, 'uint16')
    img = uint8(double(img) / 257);
elseif isa(img, 'double')
    if max(img(:)) <= 1
        img = uint8(255 * img);
    else
        img = uint8(img);
    end
elseif isa(img, 'single')
    if max(img(:)) <= 1
        img = uint8(255 * double(img));
    else
        img = uint8(img);
    end
elseif ~isa(img, 'uint8')
    img = uint8(img);
end

if ndims(img) == 2
    img = repmat(img, [1 1 3]);
elseif ndims(img) == 3
    if size(img,3) == 1
        img = repmat(img, [1 1 3]);
    elseif size(img,3) > 3
        img = img(:,:,1:3);
    end
else
    error('spt_dataset_export:InvalidImage', 'Unsupported image dimensionality.');
end
end

function save_plot(savepath)
fig = getframe(gcf);
cdata = frame2im(fig);
imwrite(cdata, savepath);
end

function close_if_valid(hFig)
if ~isempty(hFig) && isgraphics(hFig)
    close(hFig);
end
end

function delete_if_exists(fpath)
if ~isempty(fpath) && isfile(fpath)
    delete(fpath);
end
end