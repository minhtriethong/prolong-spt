function D = diffusion_coeff_v04(metadata, tracksCoordinates, opts)
% DIFFUSION_COEFF_V03_AUTOMATE (optimized)
% Rolling 2D diffusion coefficient:
%   D = MSD / (4*dt_eff)
% where dt_eff accounts for gaps between valid points.
%
% INPUTS
%   metadata.dt         : seconds per frame
%   metadata.voxel_size : physical units per pixel
%   tracksCoordinates.X : [nFrames x nTracks] (NaN for missing)
%   tracksCoordinates.Y : [nFrames x nTracks]
%
%   opts (optional)
%     .windowSize      (default 11) must be odd
%     .useParallel     (default false) parfor over tracks (can increase memory a lot on Windows)
%     .startParpool    (default false) start pool if needed
%     .storeRolling    (default true) store full rolling vectors per track (memory heavy)
%     .storeCombined   (default true) store combined non-NaN rolling values (memory heavy)
%     .outputClass     (default '') '' => single if X is single else double
%
% OUTPUT (struct D)
%   D.D_cell
%   D.D_track
%   D.D_track_rolling
%   D.D_track_rolling_combined

if nargin < 3 || isempty(opts), opts = struct(); end
opts = dflt(opts,'windowSize',11);
opts = dflt(opts,'useParallel',false);
opts = dflt(opts,'startParpool',false);
opts = dflt(opts,'storeRolling',true);
opts = dflt(opts,'storeCombined',true);
opts = dflt(opts,'outputClass','');

% ---- validate metadata ----
if ~isstruct(metadata) || ~isfield(metadata,'dt') || ~isfield(metadata,'voxel_size')
    error('metadata must contain fields: dt, voxel_size');
end
dt = double(metadata.dt);
voxel = double(metadata.voxel_size);
if ~isfinite(dt) || dt <= 0
    error('metadata.dt must be finite and > 0');
end
if ~isfinite(voxel) || voxel <= 0
    error('metadata.voxel_size must be finite and > 0');
end

% ---- validate tracksCoordinates ----
if ~isstruct(tracksCoordinates) || ~isfield(tracksCoordinates,'X') || ~isfield(tracksCoordinates,'Y')
    error('tracksCoordinates must contain fields X and Y');
end
X = tracksCoordinates.X;
Y = tracksCoordinates.Y;
if ~ismatrix(X) || ~ismatrix(Y)
    error('tracksCoordinates.X and .Y must be 2D matrices');
end
if ~isequal(size(X), size(Y))
    error('tracksCoordinates.X and .Y must have the same size');
end

[nFrames, nTracks] = size(X);

% ---- window checks ----
w = opts.windowSize;
if ~isscalar(w) || w < 3 || mod(w,2) == 0
    error('opts.windowSize must be an odd integer >= 3 (e.g., 11)');
end

% ---- output class ----
if isempty(opts.outputClass)
    if isa(X,'single') || isa(Y,'single')
        outCls = 'single';
    else
        outCls = 'double';
    end
else
    outCls = opts.outputClass;
    if ~any(strcmp(outCls, {'single','double'}))
        error('opts.outputClass must be ''single'' or ''double'' (or '''' for auto)');
    end
end

% ---- preallocate outputs ----
D_track = NaN(nTracks, 1);

storeRolling  = logical(opts.storeRolling);
storeCombined = logical(opts.storeCombined);

if storeRolling
    D_track_rolling = cell(nTracks,1);
else
    D_track_rolling = {};
end

if storeCombined && ~storeRolling
    % If you don't store rolling, we still need something to concatenate.
    D_nonNan = cell(nTracks,1);
else
    D_nonNan = {};
end

% ---- optional parallel ----
usePar = false;
if opts.useParallel
    usePar = setupPar(opts.startParpool);
end

% ---- main loop over tracks ----
if usePar
    parfor j = 1:nTracks
        Dj = diffusion_coeff_rolling_fast(X(:,j), Y(:,j), dt, voxel, w, outCls);
        D_track(j) = mean(double(Dj), 'omitnan');

        if storeRolling
            D_track_rolling{j} = Dj;
        elseif storeCombined
            v = Dj(isfinite(Dj));
            D_nonNan{j} = v(:); %#ok<PFOUS>
        end
    end
else
    for j = 1:nTracks
        Dj = diffusion_coeff_rolling_fast(X(:,j), Y(:,j), dt, voxel, w, outCls);
        D_track(j) = mean(double(Dj), 'omitnan');

        if storeRolling
            D_track_rolling{j} = Dj;
        elseif storeCombined
            v = Dj(isfinite(Dj));
            D_nonNan{j} = v(:);
        end
    end
end

% ---- aggregate outputs ----
D_cell = mean(D_track, 'omitnan');

if storeCombined
    if storeRolling
        % Preallocate combined vector to avoid giant vertcat reallocation.
        nn = cellfun(@(v) nnz(isfinite(v)), D_track_rolling);
        totalN = sum(nn);
        D_comb = NaN(totalN, 1, outCls);

        p = 1;
        for j = 1:nTracks
            v = D_track_rolling{j};
            if isempty(v), continue; end
            v = v(isfinite(v));
            m = numel(v);
            if m > 0
                D_comb(p:p+m-1) = cast(v(:), outCls);
                p = p + m;
            end
        end
    else
        % We stored only non-NaN values per track
        if isempty(D_nonNan)
            D_comb = cast([], outCls);
        else
            D_comb = vertcat(D_nonNan{:});
            D_comb = cast(D_comb, outCls);
        end
    end
else
    D_comb = cast([], outCls);
end

% ---- assemble output ----
D = struct();
D.D_cell = D_cell;
D.D_track = D_track;

if storeRolling
    D.D_track_rolling = D_track_rolling;
else
    D.D_track_rolling = {};
end

D.D_track_rolling_combined = D_comb;

end

% =================================================================
% FAST SINGLE-TRACK ROLLING D
% =================================================================
function D_rolling = diffusion_coeff_rolling_fast(X, Y, dt, voxel, windowSize, outCls)
% Computes rolling diffusion coefficient for one track using range-add + prefix sums.
% Output is length N with NaNs where undefined, matching the original behavior.

N = numel(X);
D_rolling = NaN(N,1,outCls);

half = floor(windowSize/2);
centerMin = 1 + half;
centerMax = N - half;
if centerMin > centerMax
    return; % window larger than trajectory
end

valid = isfinite(X) & isfinite(Y);
idx = find(valid);
if numel(idx) < 2
    return;
end

% positions in physical units (use double internally)
xv = double(X(idx)) * voxel;
yv = double(Y(idx)) * voxel;

dx = diff(xv);
dy = diff(yv);
d2 = dx.^2 + dy.^2;

gap = diff(idx);
dt_eff = double(gap) * dt;

goodStep = isfinite(d2) & isfinite(dt_eff) & (dt_eff > 0);
if ~any(goodStep)
    return;
end

s = idx(1:end-1);
e = idx(2:end);
s = s(goodStep);
e = e(goodStep);
Dstep = d2(goodStep) ./ (4 * dt_eff(goodStep));

% Centers i where step endpoints both lie inside [i-half, i+half]:
lo = e - half;
hi = s + half;

% Clip to valid center range (original code skips edge centers)
lo = max(lo, centerMin);
hi = min(hi, centerMax);

keep = (lo <= hi) & isfinite(Dstep);
if ~any(keep)
    return;
end

lo = lo(keep);
hi = hi(keep);
Dstep = Dstep(keep);

m = numel(Dstep);
idxAdd = [lo; hi+1];
valSum = [Dstep; -Dstep];
valCnt = [ones(m,1); -ones(m,1)];

% Range-add via difference arrays
sumDiff = accumarray(idxAdd, valSum, [N+1 1], @sum, 0);
cntDiff = accumarray(idxAdd, valCnt, [N+1 1], @sum, 0);

sumRoll = cumsum(sumDiff(1:N));
cntRoll = cumsum(cntDiff(1:N));

mask = (cntRoll > 0);
tmp = NaN(N,1);
tmp(mask) = sumRoll(mask) ./ cntRoll(mask);

% keep edges NaN to match original window-boundary behavior
tmp(1:centerMin-1) = NaN;
tmp(centerMax+1:end) = NaN;

D_rolling = cast(tmp, outCls);
end

% =================================================================
% MISC
% =================================================================
function s = dflt(s,f,v)
if ~isfield(s,f) || isempty(s.(f)), s.(f)=v; end
end

function usePar = setupPar(startPool)
usePar = false;
try
    if license('test','Distrib_Computing_Toolbox') && ~isempty(ver('parallel'))
        if isempty(gcp('nocreate'))
            if startPool
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