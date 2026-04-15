function combine_Results(Results_folderpath, nOutputSubfolders)
%COMBINE_RESULTS Usage guide for combining repeated processingInfo.mat files.
%
% Syntax
%   combine_Results(Results_folderpath, nOutputSubfolders)
%
% Purpose
%   Combine repeated processingInfo.mat files from numbered subfolders and
%   save the combined results into a new sibling output folder.
%
% Input folder format
%   Results_folderpath must point to the folder that directly contains the
%   numbered subfolders, for example:
%
%       Results/
%           subfolder-1/processingInfo.mat
%           subfolder-2/processingInfo.mat
%           ...
%
%   Important:
%   Pass the actual Results folder path, not the parent folder above it.
%
% Inputs
%   Results_folderpath
%       Character vector or scalar string containing the path to the
%       Results folder.
%
%   nOutputSubfolders
%       Number of combined output subfolders to create.
%
%       The function automatically computes:
%           repeat_no = totalNumberOfSubfolders / nOutputSubfolders
%
%       So the total number of subfolders must be divisible by
%       nOutputSubfolders.
%
% Output
%   This function saves combined files into a new sibling folder named:
%
%       <originalFolderName>_yyyymmdd_hhmmss
%
%   Example output:
%
%       Results_yyyymmdd_hhmmss/
%           subfolder-1/processingInfo.mat
%           subfolder-2/processingInfo.mat
%           ...
%
% Example
%   If the input folder contains 120 subfolders and:
%
%       nOutputSubfolders = 40
%
%   then the function will combine every 3 matching repeats:
%
%       subfolder-1  + subfolder-41 + subfolder-81  -> output subfolder-1
%       subfolder-2  + subfolder-42 + subfolder-82  -> output subfolder-2
%
% Requirement
%   The function combine_processingInfo(inputPaths, outputMatPath) must
%   already be available on the MATLAB path.
%
% Notes
%   This function does not return a variable. It writes the combined
%   processingInfo.mat files directly to disk.

Results_folderpath = iValidateResultsFolderPath(Results_folderpath);
nOutputSubfolders = iValidateNOutputSubfolders(nOutputSubfolders);

[parentDir, resultsFolderName] = fileparts(Results_folderpath);
if isempty(resultsFolderName)
    [parentDir, resultsFolderName] = fileparts(parentDir);
end

subfolders = iGetSortedSubfolders(Results_folderpath);
n = numel(subfolders);

if n == 0
    error('combine_Results:NoSubfolders', ...
        'No valid subfolders were found in: %s', Results_folderpath);
end

if nOutputSubfolders > n
    error('combine_Results:InvalidNOutputSubfolders', ...
        'nOutputSubfolders (%d) cannot exceed the number of subfolders (%d).', ...
        nOutputSubfolders, n);
end

if mod(n, nOutputSubfolders) ~= 0
    error('combine_Results:SubfolderCountMismatch', ...
        ['The number of subfolders (%d) is not divisible by nOutputSubfolders (%d). ' ...
        'Cannot form equal groups.'], n, nOutputSubfolders);
end

repeat_no = n / nOutputSubfolders;

if repeat_no < 2
    error('combine_Results:InvalidDerivedRepeatNo', ...
        ['Derived repeat_no is %d. At least 2 inputs are required for ' ...
        'combine_processingInfo.'], repeat_no);
end

timestamp = iMakeTimestamp();
outputRoot = fullfile(parentDir, [resultsFolderName '_' timestamp]);
outputRoot = iMakeUniquePath(outputRoot);

mkdir(outputRoot);

fprintf('Input folder : %s\n', Results_folderpath);
fprintf('Subfolders   : %d\n', n);
fprintf('nOutputSubfolders: %d\n', nOutputSubfolders);
fprintf('repeat_no         : %d\n', repeat_no);
fprintf('Output folder: %s\n\n', outputRoot);

subfolderNames = {subfolders.name};
matPaths = cell(n, 1);

for i = 1:n
    matPaths{i} = fullfile(Results_folderpath, subfolderNames{i}, 'processingInfo.mat');
    if ~isfile(matPaths{i})
        error('combine_Results:MissingProcessingInfo', ...
            'Missing file: %s', matPaths{i});
    end
end

for k = 1:nOutputSubfolders
    inputIdx = k + (0:repeat_no-1) * nOutputSubfolders;
    inputPaths = matPaths(inputIdx);

    outputSubfolder = fullfile(outputRoot, subfolderNames{k});
    if ~isfolder(outputSubfolder)
        mkdir(outputSubfolder);
    end

    outputMatPath = fullfile(outputSubfolder, 'processingInfo.mat');

    fprintf('[%d/%d] Combining:\n', k, nOutputSubfolders);
    for j = 1:numel(inputPaths)
        fprintf('    %s\n', inputPaths{j});
    end
    fprintf(' -> %s\n\n', outputMatPath);

    combine_processingInfo(inputPaths, outputMatPath);
end

fprintf('Done.\n');
end


function folderPath = iValidateResultsFolderPath(folderPath)
if isstring(folderPath)
    if ~isscalar(folderPath)
        error('combine_Results:InvalidFolderPath', ...
            'Results_folderpath must be a character vector or scalar string.');
    end
    folderPath = char(folderPath);
end

if ~(ischar(folderPath) && isrow(folderPath))
    error('combine_Results:InvalidFolderPath', ...
        'Results_folderpath must be a character vector or scalar string.');
end

folderPath = strtrim(folderPath);

if isempty(folderPath)
    error('combine_Results:EmptyFolderPath', ...
        'Results_folderpath cannot be empty.');
end

if ~isfolder(folderPath)
    error('combine_Results:FolderNotFound', ...
        'Folder not found: %s', folderPath);
end

[~, currentName] = fileparts(folderPath);
if isempty(currentName)
    [~, currentName] = fileparts(fileparts(folderPath));
end

candidateResultsPath = fullfile(folderPath, 'Results');
if ~strcmpi(currentName, 'Results') && isfolder(candidateResultsPath)
    listing = dir(folderPath);
    listing = listing([listing.isdir]);
    names = {listing.name};
    keep = ~ismember(names, {'.', '..'});
    names = names(keep);
    keep = cellfun(@(x) isempty(x) || x(1) ~= '.', names);
    names = names(keep);

    hasNumericSubfoldersHere = any(~cellfun(@isempty, ...
        regexp(names, '(\d+)(?!.*\d)', 'once')));

    if ~hasNumericSubfoldersHere
        warning('combine_Results:AdjustedResultsFolderPath', ...
            ['Results_folderpath looks like a parent folder, not the actual Results folder. ' ...
             'Please pass this path instead: %s. Using it automatically now.'], ...
            candidateResultsPath);
        folderPath = candidateResultsPath;
    end
end
end


function nOutputSubfolders = iValidateNOutputSubfolders(nOutputSubfolders)
if ~(isnumeric(nOutputSubfolders) && isscalar(nOutputSubfolders) ...
        && isfinite(nOutputSubfolders) && nOutputSubfolders >= 1 ...
        && mod(nOutputSubfolders, 1) == 0)
    error('combine_Results:InvalidNOutputSubfolders', ...
        'nOutputSubfolders must be an integer >= 1.');
end

nOutputSubfolders = double(nOutputSubfolders);
end


function subfolders = iGetSortedSubfolders(rootFolder)
listing = dir(rootFolder);
listing = listing([listing.isdir]);

names = {listing.name};

keep = ~ismember(names, {'.', '..'});
listing = listing(keep);
names = {listing.name};

% Ignore hidden folders such as .git, .svn, etc.
keep = cellfun(@(x) isempty(x) || x(1) ~= '.', names);
listing = listing(keep);

if isempty(listing)
    subfolders = listing;
    return;
end

numericIds = zeros(numel(listing), 1);
sortKeys = cell(numel(listing), 1);

for i = 1:numel(listing)
    numericIds(i) = iExtractLastNumber(listing(i).name);
    sortKeys{i} = sprintf('%020.0f_%s', numericIds(i), lower(listing(i).name));
end

[~, order] = sort(sortKeys);
subfolders = listing(order);
end


function n = iExtractLastNumber(folderName)
token = regexp(folderName, '(\d+)(?!.*\d)', 'tokens', 'once');

if isempty(token)
    error('combine_Results:InvalidSubfolderName', ...
        ['Subfolder name must contain a numeric suffix for ordering. ' ...
        'Problematic folder: %s'], folderName);
end

n = str2double(token{1});

if isnan(n)
    error('combine_Results:InvalidSubfolderName', ...
        'Could not parse numeric suffix from folder name: %s', folderName);
end
end


function timestamp = iMakeTimestamp()
c = clock; % [year month day hour minute seconds]
timestamp = sprintf('%04d%02d%02d_%02d%02d%02d', ...
    c(1), c(2), c(3), c(4), c(5), floor(c(6)));
end


function outPath = iMakeUniquePath(basePath)
outPath = basePath;
counter = 1;

while isfolder(outPath) || isfile(outPath)
    outPath = sprintf('%s_%d', basePath, counter);
    counter = counter + 1;
end
end