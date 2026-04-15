function processingInfo = combine_processingInfo(inputPaths, outputMatPath)
%COMBINE_PROCESSINGINFO Combine two or more processingInfo.mat files.
%
%   processingInfo = COMBINE_PROCESSINGINFO(inputPaths)
%   processingInfo = COMBINE_PROCESSINGINFO(inputPaths, outputMatPath)
%
% INPUTS
%   inputPaths    Cell array or string array of paths to processingInfo.mat
%                 files. Folder paths are also accepted; in that case this
%                 function will look for "processingInfo.mat" inside each
%                 folder.
%
%   outputMatPath Optional output .mat path. If this is a folder, the file
%                 "processingInfo_combined.mat" will be written inside it.
%
% MERGE RULES
%   1) Concatenate vertically ("below"):
%        processingInfo.Coords
%        processingInfo.track_features.D.D_track
%        processingInfo.track_features.D.D_track_rolling
%        processingInfo.track_features.D.D_track_rolling_combined
%
%   2) Concatenate horizontally ("to the right"):
%        processingInfo.tracksCoordinates.X
%        processingInfo.tracksCoordinates.Y
%        processingInfo.track_features.travelDistances
%        processingInfo.track_features.travelTime
%        processingInfo.track_features.majorAxis
%        processingInfo.track_features.minorAxis
%
%   3) Weighted average for:
%        processingInfo.track_features.D.D_cell
%      using the number of tracks in track_features.D.D_track as weights.
%
%   4) All other fields are kept from the first input file exactly as-is.
%
% EXAMPLE
%   files = {
%       'C:\data\processingInfo-1.mat'
%       'C:\data\processingInfo-2.mat'
%       'C:\data\processingInfo-3.mat'
%   };
%
%   processingInfo = combine_processingInfo(files, 'C:\data\processingInfo_combined.mat');
%
% NOTES
%   - The first input is used as the template.
%   - This function is written to combine more than two files by merging
%     them sequentially.
%   - By request, fields not listed above are NOT recalculated.

narginchk(1, 2);

if nargin < 2
    outputMatPath = '';
end

inputPaths = iNormalizeInputPaths(inputPaths);

if numel(inputPaths) < 2
    error('combine_processingInfo:NeedAtLeastTwoInputs', ...
        'Provide at least two processingInfo.mat inputs.');
end

processingInfo = [];
nValid = 0;

for k = 1:numel(inputPaths)
    try
        nextProcessingInfo = iLoadProcessingInfo(inputPaths{k});

        if nValid == 0
            processingInfo = nextProcessingInfo;
            nValid = 1;
        else
            processingInfo = iMergeTwoProcessingInfo(processingInfo, nextProcessingInfo);
            nValid = nValid + 1;
        end
    catch ME
        warning('combine_processingInfo:SkippedInput', ...
            'Skipping processingInfo file: %s. Reason: %s', ...
            inputPaths{k}, ME.message);
    end
end

if nValid == 0
    error('combine_processingInfo:NoValidInputs', ...
        'No valid processingInfo.mat files were found to combine.');
elseif nValid == 1
    warning('combine_processingInfo:OnlyOneValidInput', ...
        'Only one valid processingInfo.mat file was found. Returning it without merging.');
end

if ~isempty(outputMatPath)
    outputMatPath = char(string(outputMatPath));

    if isfolder(outputMatPath)
        outputMatPath = fullfile(outputMatPath, 'processingInfo_combined.mat');
    end

    outDir = fileparts(outputMatPath);
    if ~isempty(outDir) && ~isfolder(outDir)
        mkdir(outDir);
    end

    iSaveProcessingInfo(outputMatPath, processingInfo);
end
end


function inputPaths = iNormalizeInputPaths(inputPaths)
if ischar(inputPaths)
    inputPaths = {inputPaths};
elseif isstring(inputPaths)
    inputPaths = cellstr(inputPaths(:));
elseif iscell(inputPaths)
    inputPaths = inputPaths(:);
else
    error('combine_processingInfo:InvalidInputPaths', ...
        'inputPaths must be a char array, string array, or cell array of paths.');
end

for k = 1:numel(inputPaths)
    p = inputPaths{k};

    if ~(ischar(p) || (isstring(p) && isscalar(p)))
        error('combine_processingInfo:InvalidPathEntry', ...
            'Each entry in inputPaths must be a character vector or scalar string.');
    end

    p = char(string(p));

    if isfolder(p)
        p = fullfile(p, 'processingInfo.mat');
    end

    if ~isfile(p)
        error('combine_processingInfo:FileNotFound', ...
            'Could not find file: %s', p);
    end

    inputPaths{k} = p;
end
end


function processingInfo = iLoadProcessingInfo(matPath)
S = load(matPath, 'processingInfo');

if ~isfield(S, 'processingInfo')
    error('combine_processingInfo:MissingVariable', ...
        'File does not contain variable "processingInfo": %s', matPath);
end

processingInfo = S.processingInfo;

if ~isstruct(processingInfo)
    error('combine_processingInfo:InvalidProcessingInfo', ...
        'Variable "processingInfo" in %s must be a struct.', matPath);
end

% Required fields for the requested merge operations.
iRequireNestedField(processingInfo, {'Coords'});
iRequireNestedField(processingInfo, {'tracksCoordinates','X'});
iRequireNestedField(processingInfo, {'tracksCoordinates','Y'});
iRequireNestedField(processingInfo, {'track_features','D','D_track'});
iRequireNestedField(processingInfo, {'track_features','D','D_track_rolling'});
iRequireNestedField(processingInfo, {'track_features','D','D_track_rolling_combined'});
iRequireNestedField(processingInfo, {'track_features','D','D_cell'});
iRequireNestedField(processingInfo, {'track_features','travelDistances'});
iRequireNestedField(processingInfo, {'track_features','travelTime'});
iRequireNestedField(processingInfo, {'track_features','majorAxis'});
iRequireNestedField(processingInfo, {'track_features','minorAxis'});
end


function out = iMergeTwoProcessingInfo(a, b)
out = a;

nA = iTrackCount(a);
nB = iTrackCount(b);

% 1) Put below
out.Coords = iConcatBelow(a.Coords, b.Coords);
out.track_features.D.D_track = iConcatBelow(a.track_features.D.D_track, b.track_features.D.D_track);
out.track_features.D.D_track_rolling = iConcatBelow(a.track_features.D.D_track_rolling, b.track_features.D.D_track_rolling);
out.track_features.D.D_track_rolling_combined = iConcatBelow( ...
    a.track_features.D.D_track_rolling_combined, ...
    b.track_features.D.D_track_rolling_combined);

% 2) Put to the right
out.tracksCoordinates.X = iConcatRight(a.tracksCoordinates.X, b.tracksCoordinates.X);
out.tracksCoordinates.Y = iConcatRight(a.tracksCoordinates.Y, b.tracksCoordinates.Y);
out.track_features.travelDistances = iConcatRight(a.track_features.travelDistances, b.track_features.travelDistances);
out.track_features.travelTime = iConcatRight(a.track_features.travelTime, b.track_features.travelTime);
out.track_features.majorAxis = iConcatRight(a.track_features.majorAxis, b.track_features.majorAxis);
out.track_features.minorAxis = iConcatRight(a.track_features.minorAxis, b.track_features.minorAxis);

% 3) Weighted D_cell using D_track lengths as weights.
out.track_features.D.D_cell = iWeightedMean(a.track_features.D.D_cell, nA, ...
    b.track_features.D.D_cell, nB);

% 4) All other fields remain from the first input (already satisfied by
%    starting from out = a and only overwriting the requested fields).
end


function n = iTrackCount(processingInfo)
dTrack = processingInfo.track_features.D.D_track;

if isempty(dTrack)
    n = 0;
else
    n = length(dTrack);
end
end


function out = iConcatBelow(a, b)
if isempty(a)
    out = b;
    return;
end

if isempty(b)
    out = a;
    return;
end

if istable(a) && istable(b)
    out = iConcatTablesVertically(a, b);
    return;
end

if isvector(a) && isvector(b)
    out = [a(:); b(:)];
    return;
end

try
    out = [a; b];
catch ME
    error('combine_processingInfo:ConcatBelowFailed', ...
        'Vertical concatenation failed: %s', ME.message);
end
end


function out = iConcatRight(a, b)
if isempty(a)
    out = b;
    return;
end

if isempty(b)
    out = a;
    return;
end

if isvector(a) && isvector(b)
    out = [a(:).', b(:).'];
    return;
end

if size(a, 1) ~= size(b, 1)
    error('combine_processingInfo:ConcatRightRowMismatch', ...
        'Horizontal concatenation failed because row counts differ (%d vs %d).', ...
        size(a, 1), size(b, 1));
end

try
    out = [a, b];
catch ME
    error('combine_processingInfo:ConcatRightFailed', ...
        'Horizontal concatenation failed: %s', ME.message);
end
end


function out = iConcatTablesVertically(a, b)
varsA = a.Properties.VariableNames;
varsB = b.Properties.VariableNames;

if numel(varsA) ~= numel(varsB) || ~all(ismember(varsA, varsB))
    error('combine_processingInfo:CoordsTableMismatch', ...
        'Coords tables do not have matching variable names.');
end

% Reorder b to match a before concatenation.
b = b(:, varsA);

out = [a; b];
end


function out = iWeightedMean(d1, n1, d2, n2)
if isempty(n1), n1 = 0; end
if isempty(n2), n2 = 0; end

n1 = double(n1);
n2 = double(n2);
d1 = double(d1);
d2 = double(d2);

if n1 <= 0 && n2 <= 0
    out = d1;
elseif n1 <= 0
    out = d2;
elseif n2 <= 0
    out = d1;
else
    out = (n1 .* d1 + n2 .* d2) ./ (n1 + n2);
end
end


function iSaveProcessingInfo(outputMatPath, processingInfo)
try
    save(outputMatPath, 'processingInfo');
catch
    save(outputMatPath, 'processingInfo', '-v7.3');
end
end


function iRequireNestedField(S, pathParts)
cur = S;

for k = 1:numel(pathParts)
    thisField = pathParts{k};

    if ~isstruct(cur) || ~isfield(cur, thisField)
        error('combine_processingInfo:MissingField', ...
            'Missing required field: %s', strjoin(pathParts, '.'));
    end

    cur = cur.(thisField);
end

if isempty(cur)
    error('combine_processingInfo:EmptyField', ...
        'Required field is empty: %s', strjoin(pathParts, '.'));
end
end
