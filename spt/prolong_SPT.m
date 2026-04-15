function prolong_SPT(inputFolder, boxSize, minNetGradient, voxelSize, frameInterval, contrastMinVal, contrastMaxVal)

narginchk(5, 7);

if nargin < 6
    contrastMinVal = [];
end
if nargin < 7
    contrastMaxVal = [];
end

inputFolder = char(string(inputFolder));
boxSize = normalizeBoxSize(boxSize);
minNetGradient = double(minNetGradient);
voxelSize = double(voxelSize);
frameInterval = double(frameInterval);

validateattributes(boxSize, {'numeric'}, {'scalar','finite','real','integer','>=',1}, mfilename, 'boxSize');
if mod(double(boxSize), 2) ~= 1
    error('prolong_SPT:InvalidBoxSize', 'Box size must be an odd integer.');
end

validateattributes(minNetGradient, {'numeric'}, {'scalar','finite','real','>=',0}, mfilename, 'minNetGradient');
validateattributes(voxelSize, {'numeric'}, {'scalar','finite','real','>',0}, mfilename, 'voxelSize');
validateattributes(frameInterval, {'numeric'}, {'scalar','finite','real','>',0}, mfilename, 'frameInterval');

[contrastMinVal, contrastMaxVal] = normalizeContrast(contrastMinVal, contrastMaxVal);

% Match the GUI default checkbox value.
% Set this to false to disable Gaussian refinement.
useGaussianFit = true;

if ~isfolder(inputFolder)
    error('prolong_SPT:MissingInputFolder', 'Input folder does not exist: %s', inputFolder);
end

rootFolder = resolveRootFolder();
if ~isempty(rootFolder) && isfolder(rootFolder)
    addpath(genpath(rootFolder));
end

ensureBioFormatsInitialized(rootFolder);

requiredFcns = { ...
    'bfGetReader', ...
    'bfGetPlane', ...
    'single_particle_detection_standalone_v01', ...
    'single_particle_detection_standalone_subpixel_quadratic_v01', ...
    'particle_tracking_utrack_v03'};

if useGaussianFit
    requiredFcns{end+1} = 'single_particle_detection_standalone_subpixel_precision_v01';
end

missing = requiredFcns(cellfun(@(f) exist(f, 'file') ~= 2, requiredFcns));
if ~isempty(missing)
    error('prolong_SPT:MissingDependencies', ...
        'Missing required functions on the MATLAB path: %s', strjoin(missing, ', '));
end

results_folder = chooseResultsFolder(inputFolder);
[ok, msg] = mkdir(results_folder);
if ~ok && ~isfolder(results_folder)
    error('prolong_SPT:ResultsFolderCreateFailed', ...
        'Could not create results folder: %s\n%s', results_folder, msg);
end

logPath = fullfile(results_folder, 'NNB_log.txt');
fid = fopen(logPath, 'a');
if fid < 0
    fid = [];
end
fidCleanup = onCleanup(@() safeFclose(fid)); %#ok<NASGU>

status(sprintf('Input folder: %s', inputFolder));
status(sprintf('Results folder: %s', results_folder));

logLine(fid, sprintf('=== Run %s ===', datestr(now)));
logLine(fid, sprintf('InputFolder: %s', inputFolder));
logLine(fid, sprintf(['BoxSize=%g, MinNetGradient=%g, VoxelSize=%g, ' ...
    'FrameInterval=%g, DetectOnAllSlices=%d, UseGaussianSubpixelPrecision=%d'], ...
    boxSize, minNetGradient, voxelSize, frameInterval, 1, double(useGaussianFit)));

status('Scanning input folder for Bio-Formats-readable files...');
files = listReadableFiles(inputFolder);

totalAcquisitions = numel(files);
processedCount = 0;

for index = 1:totalAcquisitions
    currentfilename = files(index).name;
    sourceImagePath = fullfile(files(index).folder, currentfilename);
    [~, imageName, ~] = fileparts(currentfilename);

    status(sprintf('Processing %s (%d/%d)', imageName, index, totalAcquisitions));
    logLine(fid, sprintf('INFO [%s]: started (%d/%d).', imageName, index, totalAcquisitions));

    subFolderPath = fullfile(results_folder, imageName);
    if ~isfolder(subFolderPath)
        [ok, msg] = mkdir(subFolderPath);
        if ~ok
            logLine(fid, sprintf('ERROR [%s]: mkdir failed for %s: %s', imageName, subFolderPath, msg));
            continue;
        end
    end

    staleFigPath = fullfile(subFolderPath, 'interactiveFig.fig');
    try
        if isfile(staleFigPath)
            delete(staleFigPath);
        end
    catch
    end

    processingInfo = struct();
    processingInfo.AppInputSettings = struct( ...
        'InputFolder', inputFolder, ...
        'BoxSize', boxSize, ...
        'MinNetGradient', minNetGradient, ...
        'VoxelSize', voxelSize, ...
        'FrameInterval', frameInterval, ...
        'ContrastMinVal', contrastMinVal, ...
        'ContrastMaxVal', contrastMaxVal, ...
        'DetectOnAllSlices', true, ...
        'UseGaussianSubpixelPrecision', logical(useGaussianFit));
    processingInfo.imageName = imageName;
    processingInfo.sourceImagePath = sourceImagePath;
    processingInfo.status = "initializing";

    infoPath = fullfile(subFolderPath, 'processingInfo.mat');
    safeSaveMat(infoPath, processingInfo);

    try
        processingInfo.status = "reading movie";
        safeSaveMat(infoPath, processingInfo);

        reader = bfGetReader(sourceImagePath);
        readerCleanup = onCleanup(@() safeCloseReader(reader)); %#ok<NASGU>

        sizeX = reader.getSizeX();
        sizeY = reader.getSizeY();
        sizeZ = reader.getSizeZ();
        sizeC = reader.getSizeC();
        sizeT = reader.getSizeT();

        if sizeT < 1
            error('prolong_SPT:InvalidSizeT', ...
                'Bio-Formats reports sizeT < 1 for file: %s', sourceImagePath);
        end

        if sizeT > 1 && sizeZ > 1
            error('prolong_SPT:Unsupported4DStack', ...
                ['File %s appears to be a 4D stack (Z=%d and T=%d). ' ...
                'This pipeline supports only 2D images and 3D stacks ' ...
                '(either Z or T, not both).'], ...
                sourceImagePath, sizeZ, sizeT);
        end

        cIdx = 1;
        if sizeC < 1
            error('prolong_SPT:InvalidSizeC', ...
                'Bio-Formats reports sizeC < 1 for file: %s', sourceImagePath);
        end

        frameDim = 'T';
        nFrames = sizeT;
        if sizeT == 1 && sizeZ > 1
            frameDim = 'Z';
            nFrames = sizeZ;
            logLine(fid, sprintf('INFO [%s]: sizeT=1 and sizeZ=%d, treating Z as frames.', imageName, sizeZ));
        end

        planeData = zeros(sizeY, sizeX, nFrames, 'uint16');
        progressStride = max(1, round(nFrames / 10));

        for f = 1:nFrames
            if frameDim == 'T'
                zIdx = 1;
                tIdx = f;
            else
                zIdx = f;
                tIdx = 1;
            end

            iPlane = reader.getIndex(zIdx - 1, cIdx - 1, tIdx - 1) + 1;
            planeData(:, :, f) = bfGetPlane(reader, iPlane);

            if f == 1 || f == nFrames || mod(f, progressStride) == 0
                status(sprintf('Reading %s (%d/%d)', imageName, f, nFrames));
            end
        end

        processingInfo.status = "movie loaded";
        processingInfo.sizeX = sizeX;
        processingInfo.sizeY = sizeY;
        processingInfo.sizeZ = sizeZ;
        processingInfo.sizeC = sizeC;
        processingInfo.sizeT = sizeT;
        processingInfo.frameDimUsed = frameDim;
        processingInfo.nFramesUsed = nFrames;
        safeSaveMat(infoPath, processingInfo);

        if nFrames < 3
            error('prolong_SPT:TooFewFrames', ...
                'File %s has only %d frames; need at least 3.', currentfilename, nFrames);
        end

        processingInfo.status = "detecting particles";
        safeSaveMat(infoPath, processingInfo);

        status(sprintf('Detecting particles in %s', imageName));
        centers = single_particle_detection_standalone_v01(planeData, boxSize);

        if isempty(centers)
            logLine(fid, sprintf('WARNING [%s]: no centers found.', imageName));
            processingInfo.status = "no centers found";
            safeSaveMat(infoPath, processingInfo);
            continue;
        end

        if size(centers, 2) < 4
            error('prolong_SPT:InvalidCenters', ...
                'Centers output for %s has %d columns; expected >= 4.', ...
                imageName, size(centers, 2));
        end

        centers_filtered = centers(centers(:, 4) > minNetGradient, :);

        if isempty(centers_filtered)
            logLine(fid, sprintf('WARNING [%s]: all centers filtered out (minNetGradient=%g).', ...
                imageName, minNetGradient));
            processingInfo.status = "no centers after filtering";
            safeSaveMat(infoPath, processingInfo);
            continue;
        end

        status(sprintf('Fitting particles in %s (quadratic)', imageName));
        centers_subpixel_quadratic_all = ...
            single_particle_detection_standalone_subpixel_quadratic_v01(planeData, centers_filtered);

        if useGaussianFit
            status(sprintf('Refining particles in %s (Gaussian fit)', imageName));

            fitFrames = double(centers_subpixel_quadratic_all.frame);
            fitX = double(centers_subpixel_quadratic_all.x);
            fitY = double(centers_subpixel_quadratic_all.y);
            fitNetGradient = double(centers_subpixel_quadratic_all.net_gradient);

            fitRows = find(isfinite(fitFrames) & isfinite(fitX) & isfinite(fitY));
            centers_subpixel_precision_all = centers_subpixel_quadratic_all;

            if ~isempty(fitRows)
                centers_for_gaussian = [ ...
                    fitFrames(fitRows), ...
                    round(fitX(fitRows)), ...
                    round(fitY(fitRows)), ...
                    fitNetGradient(fitRows)];

                gaussian_locs = single_particle_detection_standalone_subpixel_precision_v01( ...
                    planeData, centers_for_gaussian, boxSize);

                if ~istable(gaussian_locs) || height(gaussian_locs) ~= numel(fitRows)
                    error('prolong_SPT:GaussianFitMismatch', ...
                        'Gaussian subpixel precision output size mismatch.');
                end

                goodMask = logical(gaussian_locs.success);
                if any(goodMask)
                    goodRows = fitRows(goodMask);
                    centers_subpixel_precision_all.x(goodRows) = single(gaussian_locs.x(goodMask));
                    centers_subpixel_precision_all.y(goodRows) = single(gaussian_locs.y(goodMask));

                    if ismember('success', centers_subpixel_precision_all.Properties.VariableNames)
                        centers_subpixel_precision_all.success(goodRows) = true;
                    end
                end
            end
        else
            centers_subpixel_precision_all = centers_subpixel_quadratic_all;
        end

        processingInfo.Coords = centers_subpixel_precision_all;
        processingInfo.numDetectedParticles = size(centers_subpixel_precision_all, 1);

        processingInfo.status = "tracking particles";
        safeSaveMat(infoPath, processingInfo);

        status(sprintf('Tracking particles in %s', imageName));
        [~, ~, processingInfo] = particle_tracking_utrack_v03(processingInfo);

        processingInfo.status = "complete";
        safeSaveMat(infoPath, processingInfo);

        logLine(fid, sprintf('INFO [%s]: saved %d detected particles and tracking results.', ...
            imageName, processingInfo.numDetectedParticles));

        processedCount = processedCount + 1;

        clear planeData centers centers_filtered centers_subpixel_quadratic_all ...
            centers_subpixel_precision_all

    catch ME
        processingInfo.status = "error";
        processingInfo.errorMessage = ME.message;
        processingInfo.errorIdentifier = ME.identifier;
        try
            processingInfo.errorReport = getReport(ME, 'extended', 'hyperlinks', 'off');
        catch
        end
        safeSaveMat(infoPath, processingInfo);

        logLine(fid, sprintf('ERROR [%s]: %s', imageName, ME.message));
        status(sprintf('ERROR in %s: %s', imageName, ME.message));
        continue;
    end
end

msg = sprintf('Finished. Processed %d/%d files. Results: %s', ...
    processedCount, totalAcquisitions, results_folder);

status(msg);
logLine(fid, msg);
end

function boxSize = normalizeBoxSize(boxSize)
if ~(isnumeric(boxSize) && isscalar(boxSize) && isreal(boxSize) && isfinite(boxSize))
    return;
end

boxSize = max(1, round(double(boxSize)));
if mod(boxSize, 2) == 0
    boxSize = boxSize + 1;
end
end

function [contrastMinVal, contrastMaxVal] = normalizeContrast(contrastMinVal, contrastMaxVal)
if isempty(contrastMinVal) || isempty(contrastMaxVal) || ...
        ~isscalar(contrastMinVal) || ~isscalar(contrastMaxVal) || ...
        ~isfinite(contrastMinVal) || ~isfinite(contrastMaxVal)

    contrastMinVal = [];
    contrastMaxVal = [];
    return;
end

contrastMinVal = double(contrastMinVal);
contrastMaxVal = double(contrastMaxVal);

if contrastMinVal == contrastMaxVal
    contrastMaxVal = contrastMinVal + 1;
elseif contrastMinVal > contrastMaxVal
    tmp = contrastMinVal;
    contrastMinVal = contrastMaxVal;
    contrastMaxVal = tmp;
end
end

function resultsFolder = chooseResultsFolder(inputFolder)
resultsFolderBase = fullfile(inputFolder, 'Results');

if ~isfolder(resultsFolderBase) && ~isfile(resultsFolderBase)
    resultsFolder = resultsFolderBase;
    return;
end

ts = datestr(now, 'yyyymmdd_HHMMSS');
resultsFolder = makeUniqueFolder(fullfile(inputFolder, ['Results_' ts]));
end

function folder = makeUniqueFolder(baseFolder)
folder = baseFolder;
if ~isfolder(folder) && ~isfile(folder)
    return;
end

k = 1;
while true
    folder = sprintf('%s_%02d', baseFolder, k);
    if ~isfolder(folder) && ~isfile(folder)
        return;
    end
    k = k + 1;
    if k > 999
        error('prolong_SPT:UniqueFolderFailed', ...
            'Could not create a unique folder name based on: %s', baseFolder);
    end
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
        error('prolong_SPT:NoFiles', 'No files found in input folder: %s', inputFolder);
    else
        error('prolong_SPT:NoReadableFiles', ...
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

function safeFclose(fid)
try
    if ~isempty(fid) && isnumeric(fid) && fid > 2
        fclose(fid);
    end
catch
end
end

function safeSaveMat(infoPath, processingInfo)
if isempty(infoPath)
    return;
end

try
    save(infoPath, 'processingInfo');
catch
    try
        save(infoPath, 'processingInfo', '-v7.3');
    catch
        % Do not fail the full run just because a debug save failed.
    end
end
end

function logLine(fid, line)
try
    if isempty(fid) || fid < 0
        return;
    end
    fprintf(fid, '%s\n', char(string(line)));
catch
end
end

function status(msg)
fprintf('[%s] %s\n', datestr(now, 'HH:MM:SS'), char(string(msg)));
end

function ensureBioFormatsInitialized(rootFolder)
persistent done

if isempty(done) || ~done
    if bioFormatsJavaOK()
        done = true;
        return;
    end

    setupBioFormats(rootFolder);

    if ~bioFormatsJavaOK()
        error('prolong_SPT:BioFormatsInitFailed', ...
            'Bio-Formats Java classes are still unavailable after setup.');
    end

    done = true;
end
end

function tf = bioFormatsJavaOK()
tf = false;
try
    javaObject('loci.formats.ImageReader');
    tf = true;
catch
    tf = false;
end
end

function setupBioFormats(rootFolder)
bfFolder = resolveBfFolder(rootFolder);
jarPath = fullfile(bfFolder, 'bioformats_package.jar');

if ~isfolder(bfFolder)
    error('prolong_SPT:MissingBFMatlab', 'bfmatlab folder not found: %s', bfFolder);
end

addpath(genpath(bfFolder), '-begin');

if ~usejava('jvm')
    error('prolong_SPT:NoJVM', 'MATLAB JVM is not available. Bio-Formats requires Java.');
end

if ~isfile(jarPath)
    error('prolong_SPT:MissingBFJar', 'bioformats_package.jar not found: %s', jarPath);
end

dyn = string(javaclasspath('-dynamic'));
if ~any(dyn == string(jarPath))
    javaaddpath(jarPath);
end

if exist('bfCheckJavaPath', 'file') == 2
    clear bfCheckJavaPath
    bfCheckJavaPath(true);
end

javaObject('loci.formats.ImageReader');
end

function bfFolder = resolveBfFolder(rootFolder)
bfFolder = '';

candidateFolders = {};

if nargin >= 1 && ~isempty(rootFolder)
    candidateFolders{end+1} = fullfile(rootFolder, 'bfmatlab');
    candidateFolders{end+1} = rootFolder;
end

try
    readerPath = which('bfGetReader');
    if ~isempty(readerPath)
        candidateFolders{end+1} = fileparts(readerPath);
    end
catch
end

try
    thisFolder = fileparts(mfilename('fullpath'));
    candidateFolders{end+1} = fullfile(thisFolder, 'bfmatlab');
    candidateFolders{end+1} = fullfile(fileparts(thisFolder), 'bfmatlab');
catch
end

candidateFolders{end+1} = fullfile(pwd, 'bfmatlab');
candidateFolders{end+1} = fullfile(fileparts(pwd), 'bfmatlab');

for k = 1:numel(candidateFolders)
    cand = candidateFolders{k};
    if isValidBfFolder(cand)
        bfFolder = cand;
        return;
    end
end

searchStarts = {};
if nargin >= 1 && ~isempty(rootFolder)
    searchStarts{end+1} = rootFolder;
end

try
    searchStarts{end+1} = fileparts(mfilename('fullpath'));
catch
end
searchStarts{end+1} = pwd;

for k = 1:numel(searchStarts)
    bfFolder = searchUpwardsForBfFolder(searchStarts{k});
    if ~isempty(bfFolder)
        return;
    end
end

error('prolong_SPT:MissingBFMatlab', ...
    ['Could not locate bfmatlab. Expected a folder containing ' ...
    'bioformats_package.jar somewhere near prolong_SPT.m or on the MATLAB path.']);
end

function tf = isValidBfFolder(folderPath)
tf = false;
if isempty(folderPath)
    return;
end

try
    tf = isfolder(folderPath) && isfile(fullfile(folderPath, 'bioformats_package.jar'));
catch
    tf = false;
end
end

function bfFolder = searchUpwardsForBfFolder(startFolder)
bfFolder = '';

if isempty(startFolder) || ~isfolder(startFolder)
    return;
end

current = char(string(startFolder));
while true
    cand = fullfile(current, 'bfmatlab');
    if isValidBfFolder(cand)
        bfFolder = cand;
        return;
    end

    parent = fileparts(current);
    if strcmp(parent, current)
        return;
    end
    current = parent;
end
end

function rootFolder = resolveRootFolder()
rootFolder = '';

try
    thisFolder = fileparts(mfilename('fullpath'));
catch
    thisFolder = '';
end

searchStarts = {};
if ~isempty(thisFolder)
    searchStarts{end+1} = thisFolder;
end
searchStarts{end+1} = pwd;

for k = 1:numel(searchStarts)
    rootFolder = searchUpwardsForRepoRoot(searchStarts{k});
    if ~isempty(rootFolder)
        return;
    end
end

if ~isempty(thisFolder)
    rootFolder = thisFolder;
else
    rootFolder = pwd;
end
end

function rootFolder = searchUpwardsForRepoRoot(startFolder)
rootFolder = '';

if isempty(startFolder) || ~isfolder(startFolder)
    return;
end

current = char(string(startFolder));
while true
    hasBf = isfolder(fullfile(current, 'bfmatlab'));
    hasNnb = isfolder(fullfile(current, 'nnb'));
    hasSpt = isfolder(fullfile(current, 'spt'));

    if hasBf || hasNnb || hasSpt
        rootFolder = current;
        return;
    end

    parent = fileparts(current);
    if strcmp(parent, current)
        return;
    end
    current = parent;
end
end
