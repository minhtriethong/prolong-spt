function [movieInfo, featGlobalID, info] = coords_to_movieInfo_utrack(coords, opts)
% COORDS_TO_MOVIEINFO_UTRACK  Convert detections -> u-track movieInfo (2D).
%
%   [movieInfo, featGlobalID, info] = coords_to_movieInfo_utrack(coords)
%   [movieInfo, featGlobalID, info] = coords_to_movieInfo_utrack(coords, opts)
%
% PURPOSE
%   Convert per-detection rows {frame,x,y} into u-track's per-frame movieInfo struct:
%     movieInfo(f).xCoord : [N x 2]  (x, sigma_x)   sigma defaults to 0
%     movieInfo(f).yCoord : [N x 2]  (y, sigma_y)   sigma defaults to 0
%     movieInfo(f).amp    : [N x 2]  (amp, sigma_a) amp defaults to 1, sigma defaults to 0
%     movieInfo(f).num    : scalar number of detections in frame
%     movieInfo(f).allCoord : [N x 4] [x sx y sy] (2D)
%     movieInfo(f).nnDist : [N x 1] nearest-neighbor distance (for local density options)
%
% INPUT
%   coords:
%     - numeric Nx3 (or Nx>=3): [frame x y ...]
%       frame can be 0-based or 1-based; it will be normalized to start at 1.
%     - OR table with variables: frame, x, y (optional amp)
%
%   opts (struct, optional):
%     useParallel (default true)
%       Uses parfor when a parallel pool exists. Can optionally start one.
%
%     startParpool (default true)
%       If useParallel=true and no pool exists, attempt to start parpool.
%
%     useLocalDensity (default true)
%       If true, compute per-point nnDist via knnsearch(K=2).
%       If false, sets nnDist to nnDistDefault (constant) to avoid extra work.
%
%     nnDistDefault (default 1000)
%       Constant nnDist value used when useLocalDensity=false or frame has 0/1 detections.
%
%     ampColumn (default [])
%       If numeric input has an amplitude/intensity column, set this to its index
%       (e.g., 4). Otherwise amp(:,1)=1.
%
%     nFrames (default [])
%       Force number of frames in movieInfo (>= max(frame)). If empty, uses max(frame).
%
%     driftXY (default [])
%       Optional [nFrames x 2] array [dx dy] subtracted from (x,y) after frame normalization.
%       driftXY row indexing must match normalized frames (1..nFrames).
%
%     preserveOriginalRowID (default true)
%       If true, featGlobalID{f}(i) equals the ORIGINAL row index in coords/table.
%       If false, uses 1..N after filtering.
%
% OUTPUT
%   movieInfo     : struct array (length nFrames)
%   featGlobalID  : cell array (length nFrames), mapping feature index -> global detection row ID
%   info          : struct with bookkeeping:
%                    .nFrames
%                    .origMinFrame
%                    .frameShiftApplied
%                    .nDetectionsIn
%                    .nDetectionsKept
%                    .countsPerFrame
%
% NOTES
%   - This is 2D only: DO NOT add zCoord field unless you truly have z, because u-track
%     will treat it as 3D if zCoord exists.
%
% EXAMPLE
%   coords = processingInfo.Coords(:,1:3); % [frame x y]
%   [movieInfo, featGlobalID, info] = coords_to_movieInfo_utrack(coords, struct( ...
%       'useParallel', true, ...
%       'useLocalDensity', true));
%

if nargin < 2 || isempty(opts), opts = struct(); end
opts = dflt(opts,'useParallel',false);
opts = dflt(opts,'startParpool',true);
opts = dflt(opts,'useLocalDensity',true);
opts = dflt(opts,'nnDistDefault',1000);
opts = dflt(opts,'ampColumn',[]);
opts = dflt(opts,'nFrames',[]);
opts = dflt(opts,'driftXY',[]);
opts = dflt(opts,'preserveOriginalRowID',true);

% ---- Parse input: frame,x,y (+ optional amp) and original row IDs ----
[frame, x, y, amp, origRow] = parseCoords(coords, opts.ampColumn);

info = struct();
info.nDetectionsIn = numel(frame);

% filter non-finite
good = isfinite(frame) & isfinite(x) & isfinite(y);
frame = frame(good);
x = x(good);
y = y(good);
amp = amp(good);
origRow = origRow(good);

info.nDetectionsKept = numel(frame);

if isempty(frame)
    % No detections: return empty movieInfo
    movieInfo = repmat(struct('xCoord',[],'yCoord',[],'amp',[],'num',0,'allCoord',[],'nnDist',[]), 0, 1);
    featGlobalID = {};
    info.nFrames = 0;
    info.origMinFrame = NaN;
    info.frameShiftApplied = NaN;
    info.countsPerFrame = [];
    return;
end

% frame must be integer-like
if any(abs(frame - round(frame)) > 0)
    error('Frame indices must be integers (or integer-like).');
end
frame = round(frame);

% ---- Normalize frames to start at 1 ----
origMin = min(frame);
if origMin == 0
    frame = frame + 1; % 0-based -> 1-based
end
origMin2 = min(frame); % after 0->1 conversion
shift = origMin2 - 1;
frame = frame - shift; % now min(frame)=1

nFrames = max(frame);
if ~isempty(opts.nFrames)
    if opts.nFrames < nFrames
        error('opts.nFrames (%d) < max normalized frame (%d).', opts.nFrames, nFrames);
    end
    nFrames = opts.nFrames;
end

info.nFrames = nFrames;
info.origMinFrame = origMin;
info.frameShiftApplied = shift;

% ---- Drift correction (optional) ----
if ~isempty(opts.driftXY)
    if size(opts.driftXY,2) ~= 2 || size(opts.driftXY,1) < nFrames
        error('opts.driftXY must be [nFrames x 2] with at least %d rows.', nFrames);
    end
    x = x - opts.driftXY(frame,1);
    y = y - opts.driftXY(frame,2);
end

% ---- Local density requires knnsearch ----
useLocalDensity = logical(opts.useLocalDensity);
if useLocalDensity && exist('knnsearch','file') ~= 2
    warning('knnsearch not available. Disabling useLocalDensity and using constant nnDist.');
    useLocalDensity = false;
end

% ---- Decide parallel ----
usePar = false;
if opts.useParallel
    usePar = setupParpool(opts.startParpool);
end

% ---- Deterministic ordering: sort by (frame, origRow) ----
% This makes featGlobalID stable across runs/chunks.
[~,ord] = sortrows([frame(:), origRow(:)], [1 2]);
frame = frame(ord);
x = x(ord);
y = y(ord);
amp = amp(ord);
origRow = origRow(ord);

% Global ID mapping: preserve original row index by default
if opts.preserveOriginalRowID
    globalID = origRow;
else
    globalID = (1:numel(frame))';
end

% ---- Precompute per-frame ranges ----
counts = accumarray(frame, 1, [nFrames 1], @sum, 0);
info.countsPerFrame = counts;

startIdx = cumsum([1; counts(1:end-1)]);
endIdx   = startIdx + counts - 1;

% ---- Preallocate outputs ----
tmpl = struct('xCoord',[],'yCoord',[],'amp',[],'num',0,'allCoord',[],'nnDist',[]);
movieInfo = repmat(tmpl, nFrames, 1);
featGlobalID = cell(nFrames,1);

nnDefault = double(opts.nnDistDefault);

% ---- Build movieInfo per frame ----
if usePar
    parfor f = 1:nFrames
        n = counts(f);
        if n==0, continue; end

        i1 = startIdx(f); i2 = endIdx(f);
        xf = x(i1:i2);
        yf = y(i1:i2);
        af = amp(i1:i2);
        gid = globalID(i1:i2);

        n = numel(xf);

        xCoord = [xf, zeros(n,1)];
        yCoord = [yf, zeros(n,1)];
        ampOut = [af, zeros(n,1)];
        allCoord = [xCoord, yCoord];

        if useLocalDensity
            if n==1
                nnDist = nnDefault;
            else
                [~,d] = knnsearch([xf,yf],[xf,yf],'K',2);
                nnDist = d(:,2);
            end
        else
            nnDist = nnDefault*ones(n,1);
        end

        movieInfo(f) = struct('xCoord',xCoord,'yCoord',yCoord,'amp',ampOut, ...
                              'num',n,'allCoord',allCoord,'nnDist',nnDist);
        featGlobalID{f} = gid;
    end
else
    for f = 1:nFrames
        n = counts(f);
        if n==0, continue; end

        i1 = startIdx(f); i2 = endIdx(f);
        xf = x(i1:i2);
        yf = y(i1:i2);
        af = amp(i1:i2);
        gid = globalID(i1:i2);

        n = numel(xf);

        xCoord = [xf, zeros(n,1)];
        yCoord = [yf, zeros(n,1)];
        ampOut = [af, zeros(n,1)];
        allCoord = [xCoord, yCoord];

        if useLocalDensity
            if n==1
                nnDist = nnDefault;
            else
                [~,d] = knnsearch([xf,yf],[xf,yf],'K',2);
                nnDist = d(:,2);
            end
        else
            nnDist = nnDefault*ones(n,1);
        end

        movieInfo(f) = struct('xCoord',xCoord,'yCoord',yCoord,'amp',ampOut, ...
                              'num',n,'allCoord',allCoord,'nnDist',nnDist);
        featGlobalID{f} = gid;
    end
end

end

% ====================== local helpers ======================

function s = dflt(s, f, v)
if ~isfield(s,f) || isempty(s.(f)), s.(f)=v; end
end

function usePar = setupParpool(startParpool)
usePar = false;
try
    if license('test','Distrib_Computing_Toolbox') && ~isempty(ver('parallel'))
        if isempty(gcp('nocreate'))
            if startParpool
                parpool;
            else
                usePar = false;
                return;
            end
        end
        usePar = true;
    end
catch
    usePar = false;
end
end

function [frame,x,y,amp,origRow] = parseCoords(coords, ampColumn)
% Returns numeric column vectors: frame,x,y,amp and original row indices.

if istable(coords)
    vars = coords.Properties.VariableNames;
    if ~all(ismember({'frame','x','y'}, vars))
        error('Table input must contain variables: frame, x, y');
    end
    frame = double(coords.frame(:));
    x = double(coords.x(:));
    y = double(coords.y(:));

    if ismember('amp', vars)
        amp = double(coords.amp(:));
    else
        amp = ones(size(frame));
    end

    origRow = (1:numel(frame))';

else
    if ~isnumeric(coords) || size(coords,2) < 3
        error('Numeric input must be Nx3 or larger: [frame x y ...]');
    end
    frame = double(coords(:,1));
    x = double(coords(:,2));
    y = double(coords(:,3));

    if ~isempty(ampColumn)
        if ampColumn < 1 || ampColumn > size(coords,2)
            error('opts.ampColumn=%d is out of range for input with %d columns.', ampColumn, size(coords,2));
        end
        amp = double(coords(:,ampColumn));
    else
        amp = ones(size(frame));
    end

    origRow = (1:numel(frame))';
end

end