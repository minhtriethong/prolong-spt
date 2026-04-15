function [total_int_list, maskFolder, total_int_bg_corrected_list, area_list, intensity_density_list] = segmentation_int_measure_dataset(folderPath)
%SEGMENTATION_INT_MEASURE_DATASET
% [total_int_list, maskFolder, total_int_bg_corrected_list, area_list, intensity_density_list] = ...
%     segmentation_int_measure_dataset(folderPath)
%
% Process all Bio-Formats-readable files in folderPath:
%   1) Read only the first plane of each file (Z=1, C=1, T=1)
%   2) Huang auto-threshold that plane
%   3) Save binary mask TIFF into <folderPath>/mask
%   4) Compute total integrated intensity inside the mask
%   5) Compute background-corrected integrated intensity
%   6) Compute mask area in pixels
%   7) Read voxel size from extractND2MetadataToStruct using:
%        voxelSize = str2double(meta.file.globalMetadata.byField.dCalibration);
%   8) Compute intensity density:
%        intensityDensity = total_int_bg_corrected / (fgArea * (voxelSize^2))
%   9) Save results into total_int_list.mat and total_int_list.csv
%
% Notes:
% - Uses first channel only.
% - Only the first slice/plane is processed, even for stacks.
% - total_int_list(i), total_int_bg_corrected_list(i), area_list(i), and
%   intensity_density_list(i) match the i-th readable file in alphabetical order.
% - Existing 4-output calls remain valid; intensity_density_list is an
%   additional 5th output.

folderPath = char(string(folderPath));

if ~isfolder(folderPath)
    error('segmentation_int_measure_dataset:MissingFolder', ...
        'Input folder does not exist: %s', folderPath);
end

if exist('bfGetReader', 'file') ~= 2 || exist('bfGetPlane', 'file') ~= 2
    error('segmentation_int_measure_dataset:MissingBioFormats', ...
        'bfGetReader/bfGetPlane not found. Add bfmatlab to the MATLAB path first.');
end

if exist('extractND2MetadataToStruct', 'file') ~= 2
    error('segmentation_int_measure_dataset:MissingMetadataFunction', ...
        'extractND2MetadataToStruct not found. Add it to the MATLAB path first.');
end

maskFolder = fullfile(folderPath, 'mask');
if ~isfolder(maskFolder)
    [ok, msg] = mkdir(maskFolder);
    if ~ok
        error('segmentation_int_measure_dataset:MaskFolderCreateFailed', ...
            'Could not create mask folder: %s\n%s', maskFolder, msg);
    end
end

files = listReadableFiles(folderPath);
total_int_list = nan(numel(files), 1);
total_int_bg_corrected_list = nan(numel(files), 1);
area_list = nan(numel(files), 1);
intensity_density_list = nan(numel(files), 1);

for k = 1:numel(files)
    srcPath = fullfile(files(k).folder, files(k).name);
    [~, baseName, ~] = fileparts(files(k).name);
    outMaskPath = fullfile(maskFolder, [baseName '_mask.tif']);
    reader = [];

    try
        fprintf('[%s] Processing %s (%d/%d)\n', ...
            datestr(now, 'HH:MM:SS'), files(k).name, k, numel(files));

        % Read voxel size from metadata struct
        voxelSize = NaN;
        try
            meta = extractND2MetadataToStruct(srcPath);
            voxelSize = meta.series.physical.sizeX.value;
            if ~(isfinite(voxelSize) && voxelSize > 0)
                voxelSize = NaN;
            end
        catch ME_meta
            warning('segmentation_int_measure_dataset:VoxelSizeUnavailable', ...
                'Could not extract dCalibration from %s. intensity density will be NaN. %s', ...
                files(k).name, ME_meta.message);
        end

        reader = bfGetReader(srcPath);

        sizeZ = double(reader.getSizeZ());
        sizeC = double(reader.getSizeC());
        sizeT = double(reader.getSizeT());

        if sizeZ < 1 || sizeC < 1 || sizeT < 1
            error('Invalid Bio-Formats dimensions for file: %s', srcPath);
        end

        % Read only the first plane: Z=1, C=1, T=1
        zIdx = 1;
        cIdx = 1;
        tIdx = 1;
        iPlane = reader.getIndex(zIdx - 1, cIdx - 1, tIdx - 1) + 1;
        img = bfGetPlane(reader, iPlane);

        % Huang threshold on first plane only
        thresh = huang_threshold_value(img);
        mask = img > thresh;
        bgMask = ~mask;

        % Foreground area in pixels
        fgArea = nnz(mask);
        area_list(k) = fgArea;

        % Raw integrated intensity inside mask
        rawInt = sum(double(img(mask)));
        total_int_list(k) = rawInt;

        % Background-corrected integrated intensity
        bgArea = nnz(bgMask);
        if bgArea > 0
            bgMean = mean(double(img(bgMask)));
            scaledBg = bgMean * fgArea;
            total_int_bg_corrected_list(k) = rawInt - scaledBg;
        else
            total_int_bg_corrected_list(k) = NaN;
        end

        % Intensity density: bg-corrected intensity / foreground area in um^2
        if isfinite(total_int_bg_corrected_list(k)) && isfinite(voxelSize) && voxelSize > 0 && fgArea > 0
            intensity_density_list(k) = total_int_bg_corrected_list(k) / (fgArea * (voxelSize^2));
        else
            intensity_density_list(k) = NaN;
        end

        % Save binary mask as single-page TIFF
        write_mask_tiff(mask, outMaskPath);

        safeCloseReader(reader);
        clear img mask bgMask meta

    catch ME
        safeCloseReader(reader);

        try
            if isfile(outMaskPath)
                delete(outMaskPath);
            end
        catch
        end

        warning('segmentation_int_measure_dataset:FileFailed', ...
            'Failed on %s: %s', files(k).name, ME.message);

        total_int_list(k) = NaN;
        total_int_bg_corrected_list(k) = NaN;
        area_list(k) = NaN;
        intensity_density_list(k) = NaN;
        continue;
    end
end

save(fullfile(maskFolder, 'total_int_list.mat'), ...
    'total_int_list', 'total_int_bg_corrected_list', 'area_list', 'intensity_density_list');

write_total_int_csv( ...
    fullfile(maskFolder, 'total_int_list.csv'), ...
    total_int_list, total_int_bg_corrected_list, area_list, intensity_density_list);

fprintf('[%s] Done. Outputs saved in: %s\n', datestr(now, 'HH:MM:SS'), maskFolder);

end

function thresh = huang_threshold_value(I)
% Huang threshold on image I.
% Exact histogram is used for integer data up to 4096 gray levels.
% Wider dynamic range falls back to a 512-bin histogram.

pixels = double(I(:));
pixels = pixels(isfinite(pixels));

if isempty(pixels)
    thresh = 0;
    return;
end

minVal = min(pixels);
maxVal = max(pixels);

if minVal == maxVal
    thresh = maxVal;
    return;
end

useExact = isinteger(I) && ((maxVal - minVal) <= 4095);

if useExact
    minI = floor(minVal);
    maxI = ceil(maxVal);
    levels = minI:maxI;
    idx = round(pixels) - minI + 1;
    counts = accumarray(idx, 1, [numel(levels), 1]).';
else
    nBins = 512;
    edges = linspace(minVal, maxVal, nBins + 1);
    counts = histcounts(pixels, edges);
    levels = (edges(1:end-1) + edges(2:end)) / 2;
end

thresh = huang_from_hist(counts, levels);
end

function thresh = huang_from_hist(counts, levels)
counts = double(counts(:)).';
levels = double(levels(:)).';

nz = find(counts > 0);
if isempty(nz)
    thresh = 0;
    return;
end

first = nz(1);
last = nz(end);

if first == last
    thresh = levels(first);
    return;
end

term = 1 / max(levels(last) - levels(first), eps);

cumCount = cumsum(counts);
cumSum = cumsum(counts .* levels);

mu0 = zeros(size(counts));
mu0(first:last) = cumSum(first:last) ./ max(cumCount(first:last), eps);

tailCount = fliplr(cumsum(fliplr(counts)));
tailSum = fliplr(cumsum(fliplr(counts .* levels)));

mu1 = zeros(size(counts));
mu1(first:last-1) = tailSum(first+1:last) ./ max(tailCount(first+1:last), eps);
mu1(last) = levels(last);

bestEntropy = inf;
bestIdx = first;

for t = first:last
    entropyVal = 0;

    lowerLevels = levels(first:t);
    mu = 1 ./ (1 + term * abs(lowerLevels - mu0(t)));
    entropyVal = entropyVal + sum(counts(first:t) .* fuzzy_entropy(mu));

    if t < last
        upperLevels = levels(t+1:last);
        mu = 1 ./ (1 + term * abs(upperLevels - mu1(t)));
        entropyVal = entropyVal + sum(counts(t+1:last) .* fuzzy_entropy(mu));
    end

    if entropyVal < bestEntropy
        bestEntropy = entropyVal;
        bestIdx = t;
    end
end

thresh = levels(bestIdx);
end

function e = fuzzy_entropy(mu)
e = zeros(size(mu));
valid = (mu > 1e-6) & (mu < 0.999999);
m = mu(valid);
e(valid) = -m .* log(m) - (1 - m) .* log(1 - m);
end

function write_mask_tiff(mask, outPath)
if isfile(outPath)
    delete(outPath);
end

frame = uint8(mask);
frame(frame ~= 0) = 255;
imwrite(frame, outPath, 'tif', 'Compression', 'none');
end

function write_total_int_csv(csvPath, total_int_list, total_int_bg_corrected_list, area_list, intensity_density_list)
fid = fopen(csvPath, 'w');
if fid < 0
    error('segmentation_int_measure_dataset:CsvOpenFailed', ...
        'Could not open CSV for writing: %s', csvPath);
end

cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, 'intergrate intensity,intergrated intensity (bg corrected),area (pixel),intensity density (int/um^2)\n');

for i = 1:numel(total_int_list)
    fprintf(fid, '%.15g,%.15g,%.15g,%.15g\n', ...
        total_int_list(i), total_int_bg_corrected_list(i), area_list(i), intensity_density_list(i));
end
end

function files = listReadableFiles(inputFolder)
allFiles = dir(inputFolder);
allFiles = allFiles(~[allFiles.isdir]);

names = string({allFiles.name});
lowerNames = lower(names);

isHidden = startsWith(names, ".");
isTmp = endsWith(lowerNames, ".tmp") | endsWith(lowerNames, "~");

candidates = allFiles(~isHidden & ~isTmp);

readable = false(size(candidates));
firstME = [];
firstBadFile = '';

for k = 1:numel(candidates)
    fpath = fullfile(candidates(k).folder, candidates(k).name);
    r = [];

    try
        r = bfGetReader(fpath);
        readable(k) = true;
    catch ME
        readable(k) = false;

        if isempty(firstME)
            firstME = ME;
            firstBadFile = fpath;
        end
    end

    safeCloseReader(r);
end

files = candidates(readable);

if isempty(files)
    if isempty(firstME)
        error('segmentation_int_measure_dataset:NoFiles', ...
            'No files found in input folder: %s', inputFolder);
    else
        error('segmentation_int_measure_dataset:NoReadableFiles', ...
            ['Bio-Formats could not open any files in:\n  %s\n\n' ...
            'First failure:\n  File: %s\n  Identifier: %s\n  Message: %s'], ...
            inputFolder, firstBadFile, firstME.identifier, firstME.message);
    end
end

[~, sidx] = sort({files.name});
files = files(sidx);
end

function safeCloseReader(reader)
try
    if ~isempty(reader)
        reader.close();
    end
catch
end
end