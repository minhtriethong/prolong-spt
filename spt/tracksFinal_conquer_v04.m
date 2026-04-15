function tracksCoordinates = tracksFinal_conquer_v04(tracksFinal_divide, framelist, junction_rows, search_radius)
% tracksFinal_conquer_v04_fast
% Faster version of tracksFinal_conquer_v03.
% Key optimizations:
%   1) No dense cost(nL,nR) matrix: builds candidate edges via junction tracks only.
%   2) side_distances vectorized over B tracks (kills inner for-loop).
%   3) STORE_INTENSITY default = false (you said you don’t need it).
%
% Output format matches your expectation:
%   tracksCoordinates.X : [total_frames x nTracks]
%   tracksCoordinates.Y : [total_frames x nTracks]
%   tracksCoordinates.I : (omitted by default)

if nargin < 4 || isempty(search_radius), search_radius = 5; end

% ===================== MEMORY / PERF KNOBS =====================
USE_SINGLE      = true;      % store X/Y as single
GROW_BY_COLS    = 5000;      % grow in bigger blocks (fewer reallocs)
STORE_INTENSITY = false;     % YOU SAID YOU DON'T NEED IT
VERBOSE         = false;
MIN_FRAMES_IN_SUBTRACK_TO_KEEP = 1;  % increase to 2-5 to reduce junk tracks and speed up
% ===============================================================

cls = 'double';
if USE_SINGLE, cls = 'single'; end

% ---- basic checks ----
if ~iscell(tracksFinal_divide)
    error('tracksFinal_divide must be a cell array.');
end
if size(framelist,2) < 2
    error('framelist must be Nx2 (startFrame,endFrame).');
end

allRows  = 1:size(framelist,1);
junction_rows = unique(junction_rows(:)'); % ensure unique row indices
baseRows = setdiff(allRows, junction_rows);   % real subtracks
baseRows = sort(baseRows);
junction_rows = sort(junction_rows);

nBase = numel(baseRows);
if nBase == 0
    tracksCoordinates.X = NaN(0,0,cls);
    tracksCoordinates.Y = NaN(0,0,cls);
    if STORE_INTENSITY
        tracksCoordinates.I = NaN(0,0,cls);
    end
    return
end

total_frames = max(framelist(baseRows,2));

% Expect junction_rows roughly = nBase-1 (one bridge between each base segment pair)
nJunc = numel(junction_rows);
nLinks = min(nBase-1, nJunc);
if nJunc ~= (nBase-1) && VERBOSE
    fprintf('[conquer] warning: nBase=%d but junction_rows=%d (expected %d). Will link %d boundaries.\n', ...
        nBase, nJunc, nBase-1, nLinks);
end

% ---- Load tracks for each base subtrack ----
Xseg = cell(nBase,1); Yseg = cell(nBase,1);
Tseg = cell(nBase,1); Nseg = zeros(nBase,1);

for s = 1:nBase
    k   = baseRows(s);
    t0  = framelist(k,1);
    t1  = framelist(k,2);
    Lk  = t1 - t0 + 1;
    Tseg{s} = t0:t1;

    if k <= numel(tracksFinal_divide) && ~isempty(tracksFinal_divide{k})
        % your existing extractor
        [~,~,x,y,~] = getTracks3D_long(tracksFinal_divide{k}, Lk);
    else
        x=[]; y=[];
    end

    % Pad to Lk rows if needed
    if ~isempty(x) && size(x,1) < Lk
        pad = NaN(Lk - size(x,1), size(x,2));
        x   = [x; pad];
        y   = [y; pad];
    end

    % Optional: drop very short tracks inside each subtrack
    if MIN_FRAMES_IN_SUBTRACK_TO_KEEP > 1 && ~isempty(x)
        len = sum(~isnan(x) & ~isnan(y), 1);
        keep = (len >= MIN_FRAMES_IN_SUBTRACK_TO_KEEP);
        x = x(:,keep);
        y = y(:,keep);
    end

    if USE_SINGLE && ~isempty(x)
        x = single(x); y = single(y);
    end

    Xseg{s} = x;
    Yseg{s} = y;
    Nseg(s) = size(x,2);
end

% ---- Load junction tracks (only the ones we will actually use) ----
Xjunc = cell(nLinks,1); Yjunc = cell(nLinks,1);
Tjunc = cell(nLinks,1);
boundaryFrame = zeros(nLinks,1);

for j = 1:nLinks
    kj = junction_rows(j);
    t0 = framelist(kj,1);
    t1 = framelist(kj,2);
    Lj = t1 - t0 + 1;
    Tjunc{j} = t0:t1;

    if kj <= numel(tracksFinal_divide) && ~isempty(tracksFinal_divide{kj})
        [~,~,xj,yj,~] = getTracks3D_long(tracksFinal_divide{kj}, Lj);
    else
        xj=[]; yj=[];
    end

    if ~isempty(xj) && size(xj,1) < Lj
        pad = NaN(Lj - size(xj,1), size(xj,2));
        xj  = [xj; pad];
        yj  = [yj; pad];
    end

    if USE_SINGLE && ~isempty(xj)
        xj = single(xj); yj = single(yj);
    end

    Xjunc{j} = xj;
    Yjunc{j} = yj;

    % junction j bridges baseRows(j) -> baseRows(j+1)
    kLeft = baseRows(j);
    boundaryFrame(j) = framelist(kLeft,2);
end

% ---- Global matrices (block preallocation) ----
seg2global = cell(nBase,1);

nCols = Nseg(1);
cap   = max([nCols, GROW_BY_COLS]);

X_total = NaN(total_frames, cap, cls);
Y_total = NaN(total_frames, cap, cls);

% Seed with first segment
if nCols > 0
    X_total(Tseg{1}, 1:nCols) = Xseg{1};
    Y_total(Tseg{1}, 1:nCols) = Yseg{1};
    seg2global{1} = 1:nCols;
else
    seg2global{1} = zeros(1,0);
end

if VERBOSE
    bytesPerEl = 8; if USE_SINGLE, bytesPerEl = 4; end
    maxPossible = sum(Nseg);
    approxGB = (double(total_frames) * double(maxPossible) * bytesPerEl * 2) / 1024^3;
    fprintf('[conquer] total_frames=%d, nBase=%d, worstCaseTracks=%d, worstCase~%.1f GB (X+Y)\n', ...
        total_frames, nBase, maxPossible, approxGB);
end

% ---- Link across boundaries ----
for s = 2:nBase
    sL = s-1; sR = s;
    j  = s-1;

    % If we have no junction for this boundary, just append as new tracks
    if j > nLinks || isempty(Xjunc{j}) || isempty(Yjunc{j})
        mapR = (nCols + (1:Nseg(sR)));
        nNew = numel(mapR);

        if (nCols + nNew) > cap
            newCap = max(cap + max(GROW_BY_COLS, nNew), nCols + nNew);
            X_total(:, cap+1:newCap) = NaN;
            Y_total(:, cap+1:newCap) = NaN;
            cap = newCap;
        end

        X_total(Tseg{sR}, mapR) = Xseg{sR};
        Y_total(Tseg{sR}, mapR) = Yseg{sR};
        nCols = nCols + nNew;
        seg2global{sR} = mapR;

        if VERBOSE
            fprintf('[conquer] seg %d: no junction -> appended %d new tracks. totalCols=%d\n', ...
                s, nNew, nCols);
        end
        continue;
    end

    pairs = link_via_junction_fast( ...
        Xseg{sL}, Yseg{sL}, Tseg{sL}, ...
        Xseg{sR}, Yseg{sR}, Tseg{sR}, ...
        Xjunc{j}, Yjunc{j}, Tjunc{j}, boundaryFrame(j), search_radius);

    mapR = nan(1, Nseg(sR));
    if ~isempty(pairs)
        % vectorized mapping from right->global via left->global
        mapR(pairs(:,2)) = seg2global{sL}(pairs(:,1));
    end

    newR = find(isnan(mapR));
    nNew = numel(newR);

    if nNew > 0
        if (nCols + nNew) > cap
            newCap = max(cap + max(GROW_BY_COLS, nNew), nCols + nNew);
            X_total(:, cap+1:newCap) = NaN;
            Y_total(:, cap+1:newCap) = NaN;
            cap = newCap;
        end

        mapR(newR) = nCols + (1:nNew);
        nCols = nCols + nNew;
    end

    if ~isempty(mapR)
        X_total(Tseg{sR}, mapR) = Xseg{sR};
        Y_total(Tseg{sR}, mapR) = Yseg{sR};
    end
    seg2global{sR} = mapR;

    if VERBOSE
        fprintf('[conquer] seg %d: L=%d R=%d pairs=%d new=%d totalCols=%d\n', ...
            s, Nseg(sL), Nseg(sR), size(pairs,1), nNew, nCols);
    end
end

% Trim unused capacity
X_total = X_total(:, 1:nCols);
Y_total = Y_total(:, 1:nCols);

tracksCoordinates.X = X_total;
tracksCoordinates.Y = Y_total;

end

% =================================================================
% FAST LINKING HELPERS
% =================================================================

function pairs = link_via_junction_fast(XL,YL,TL, XR,YR,TR, XJ,YJ,TJ, b, search_radius)
% link_via_junction_fast
% Faster greedy linking:
%   - compute dL (nL x nJ), dR (nR x nJ)
%   - for each junction track k, generate only candidate (iL,iR) edges that pass radius
%   - greedy assign from sorted candidate edge list
%
% Output pairs: [iLeft jRight]

minFramesSide = 2;

nL = size(XL,2); nR = size(XR,2); nJ = size(XJ,2);
pairs = zeros(0,2);
if nL==0 || nR==0 || nJ==0, return; end

% Split junction window
leftWin  = TJ(TJ <= b);
rightWin = TJ(TJ >  b);

% Compute distances to junction tracks (vectorized over junction tracks)
dL = side_distances_fast(XL,YL,TL, XJ,YJ,TJ, leftWin,  minFramesSide);
dR = side_distances_fast(XR,YR,TR, XJ,YJ,TJ, rightWin, minFramesSide);

% Build candidate edges via junction tracks only
edgeL = [];
edgeR = [];
edgeC = [];

for k = 1:nJ
    idxL = find(dL(:,k) < search_radius);
    if isempty(idxL), continue; end
    idxR = find(dR(:,k) < search_radius);
    if isempty(idxR), continue; end

    % costs for all combinations (idxL x idxR)
    cMat = dL(idxL,k) + (dR(idxR,k))';    % [m x n]
    [LL, RR] = ndgrid(idxL, idxR);

    edgeL = [edgeL; LL(:)]; %#ok<AGROW>
    edgeR = [edgeR; RR(:)]; %#ok<AGROW>
    edgeC = [edgeC; cMat(:)]; %#ok<AGROW>
end

if isempty(edgeC)
    return;
end

% Collapse duplicate (L,R) edges by min cost (multiple junctions can propose same pair)
pairID = edgeL + (edgeR-1)*nL;
[uID, ~, ic] = unique(pairID);
cMin = accumarray(ic, edgeC, [], @min);

% Sort edges by cost
[~, order] = sort(cMin, 'ascend');
uID = uID(order);

% Greedy 1-1 assignment without building dense cost matrix
usedL = false(nL,1);
usedR = false(nR,1);

pairs = zeros(0,2);
pairsCap = min(numel(uID), min(nL,nR));
pairs(pairsCap,2) = 0; % preallocate
nPairs = 0;

for t = 1:numel(uID)
    pid = uID(t);
    iL  = mod(pid-1, nL) + 1;
    jR  = floor((pid-1)/nL) + 1;

    if ~usedL(iL) && ~usedR(jR)
        nPairs = nPairs + 1;
        pairs(nPairs,:) = [iL jR];
        usedL(iL) = true;
        usedR(jR) = true;
        if nPairs == pairsCap
            break;
        end
    end
end

pairs = pairs(1:nPairs,:);
end

function d = side_distances_fast(XA,YA,TA, XB,YB,TB, tWin, minFrames)
% side_distances_fast
% Median distances between each track in A and each track in B over tWin.
% Optimized: vectorized over B tracks (no inner for-loop over j).

nA = size(XA,2); nB = size(XB,2);
d  = inf(nA, nB);
if isempty(tWin) || nA==0 || nB==0, return; end

% TA/TB are contiguous frame vectors in your pipeline
TA0 = TA(1); TA1 = TA(end);
TB0 = TB(1); TB1 = TB(end);

% Common frames in both tracks
tCommon = tWin(tWin>=TA0 & tWin<=TA1 & tWin>=TB0 & tWin<=TB1);
if numel(tCommon) < minFrames
    return;
end

rA = tCommon - TA0 + 1;
rB = tCommon - TB0 + 1;

XA = XA(rA,:);  YA = YA(rA,:);
XB = XB(rB,:);  YB = YB(rB,:);

% For each A track, compute distances to all B tracks at once
for i = 1:nA
    xAi = XA(:,i);
    yAi = YA(:,i);

    validAi = ~(isnan(xAi) | isnan(yAi));
    if nnz(validAi) < minFrames
        continue;
    end

    % distances: [nFvalid x nB]
    dx = XB(validAi,:) - xAi(validAi);
    dy = YB(validAi,:) - yAi(validAi);
    dist = hypot(dx, dy);

    % median ignoring NaNs
    med = median(dist, 1, 'omitnan');
    cnt = sum(~isnan(dist), 1);
    med(cnt < minFrames) = inf;

    d(i,:) = med;
end
end