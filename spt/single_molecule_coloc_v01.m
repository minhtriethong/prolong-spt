function varargout = single_molecule_coloc_v01(r, coord1, coord2, coord3, varargin)
% SINGLE_MOLECULE_COLOC_V01  Mutual-nearest-neighbor colocalization for 2 or 3 coordinate lists.
%
% LEGACY USAGE (2-channel)
%   [idx1, idx2]         = single_molecule_coloc_v01(r, coord1, coord2)
%   [idx1, idx2, dist12] = single_molecule_coloc_v01(r, coord1, coord2)
%
% LEGACY USAGE (3-channel)
%   [idx1, idx2, idx3]          = single_molecule_coloc_v01(r, coord1, coord2, coord3)
%   [idx1, idx2, idx3, dist123] = single_molecule_coloc_v01(r, coord1, coord2, coord3)
%
% NAME-VALUE (NEW, backward compatible)
%   'ReturnOut' (default false)
%     When true AND you request enough outputs, an additional output struct "out"
%     is returned containing mutually exclusive bucket indices.
%
% NEW USAGE (2-channel with out)
%   [idx1, idx2, dist12, out] = single_molecule_coloc_v01(r, coord1, coord2, 'ReturnOut', true)
%
% NEW USAGE (3-channel with out)
%   [idx1, idx2, idx3, dist123, out] = single_molecule_coloc_v01(r, coord1, coord2, coord3, 'ReturnOut', true)
%
% INPUTS
%   r      : scalar radius (pixels), must be finite and >= 0
%   coord1 : [N1 x >=2] array, columns 1-2 are [X Y]
%   coord2 : [N2 x >=2] array, columns 1-2 are [X Y]
%   coord3 : optional [N3 x >=2] array, columns 1-2 are [X Y]
%
% OUTPUTS (LEGACY behavior preserved when ReturnOut=false)
%   2-image mode:
%       idx1   : indices into coord1 (original rows)
%       idx2   : indices into coord2 (original rows)
%       dist12 : distances between matched pairs (pixels)
%
%       NOTE (legacy quirk preserved):
%         If caller requests 4 outputs in 2-channel mode with ReturnOut=false,
%         this function returns: {idx1, idx2, [], dist12}
%
%   3-image mode:
%       idx1    : indices into coord1 (original rows) for strict tri-coloc
%       idx2    : indices into coord2 (original rows) for strict tri-coloc
%       idx3    : indices into coord3 (original rows) for strict tri-coloc
%       dist123 : [K x 3] distances [d12 d13 d23] for strict tri-coloc (all <= r)
%
% OUTPUT STRUCT "out" (only when ReturnOut=true AND requested)
%   2-ch:
%     out.idx12_1, out.idx12_2, out.dist12
%     out.idx1_only, out.idx2_only
%     out.counts: n12, n1_only, n2_only
%
%   3-ch (NO ambiguous bucket):
%     out.idx123_1, out.idx123_2, out.idx123_3, out.dist123
%     out.idx12_only_1, out.idx12_only_2, out.dist12_only
%     out.idx13_only_1, out.idx13_only_3, out.dist13_only
%     out.idx23_only_2, out.idx23_only_3, out.dist23_only
%     out.idx1_only, out.idx2_only, out.idx3_only
%     out.counts: n123, n12_only, n13_only, n23_only, n1_only, n2_only, n3_only
%     out.excluded: n1, n2, n3 (multi-linked conflict points not assigned to any bucket)
%
% NOTES
%   - NaN/Inf rows are ignored (based on X,Y only).
%   - Matching is one-to-one via mutual nearest neighbors within radius.
%   - Strict tri-coloc:
%       (1) coord1<->coord2 mutual NN within r
%       (2) coord1<->coord3 mutual NN within r
%       (3) coord2 and coord3 for the same coord1 also within r (d23 <= r)
%   - Pair-only buckets ("12_only", etc) require NO mutual match to the third channel.

if nargin < 3
    error('single_molecule_coloc_v01 requires at least r, coord1, coord2.');
end

% -----------------------------
% Handle optional coord3 + Name-Value pairs
% -----------------------------
nv = {};
if nargin < 4
    coord3 = [];
else
    % If the 4th input is a string/char, treat it as the start of Name-Value args
    if (ischar(coord3) || (isstring(coord3) && isscalar(coord3)))
        nv = [{coord3} varargin];
        coord3 = [];
    else
        nv = varargin;
    end
end

% Parse Name-Value
ReturnOut = false;
if ~isempty(nv)
    if mod(numel(nv), 2) ~= 0
        error('Name-value arguments must be provided as pairs.');
    end
    for k = 1:2:numel(nv)
        name = nv{k};
        val  = nv{k+1};

        if ~(ischar(name) || (isstring(name) && isscalar(name)))
            error('Name-value parameter names must be char or scalar string.');
        end

        key = lower(char(string(name)));
        switch key
            case 'returnout'
                ReturnOut = logical(val);
                if ~isscalar(ReturnOut)
                    error('ReturnOut must be a scalar logical.');
                end
            otherwise
                error('Unknown name-value parameter: %s', char(string(name)));
        end
    end
end

% Validate r
if ~isscalar(r) || ~isnumeric(r) || ~isfinite(r) || r < 0
    error('r must be a finite, non-negative numeric scalar.');
end
r2 = r * r;

% Clean coordinates (take first 2 columns, drop non-finite rows)
[c1, orig1] = cleanCoords(coord1);
[c2, orig2] = cleanCoords(coord2);

use3 = ~isempty(coord3);
if use3
    [c3, orig3] = cleanCoords(coord3);
end

% Decide whether we will actually build out
wantOut = false;
if ~use3
    wantOut = ReturnOut && (nargout >= 4);
else
    wantOut = ReturnOut && (nargout >= 5);
end

% =========================
% 2-channel mode
% =========================
if ~use3
    [i1v, i2v, d12] = mutualNN_withinRadius(c1, c2, r2);

    idx1 = orig1(i1v);
    idx2 = orig2(i2v);

    % ----- Legacy outputs (ReturnOut=false OR out not requested) -----
    if ~wantOut
        if nargout <= 2
            varargout = {idx1, idx2};
        elseif nargout == 3
            varargout = {idx1, idx2, d12};
        else
            % legacy bug-proof: some callers ask 4 outputs; keep dist in 4th
            varargout = {idx1, idx2, [], d12};
        end
        % pad for callers asking too many outputs
        if numel(varargout) < nargout
            varargout(end+1:nargout) = {[]};
        end
        return;
    end

    % ----- Extended output struct (ReturnOut=true AND requested) -----
    out = struct();
    out.mode = '2ch';
    out.r = r;

    out.idx12_1 = idx1(:);
    out.idx12_2 = idx2(:);
    out.dist12  = d12(:);

    n1v = size(c1,1);
    n2v = size(c2,1);

    used1 = false(n1v,1); used1(i1v) = true;
    used2 = false(n2v,1); used2(i2v) = true;

    i1_only_v = find(~used1);
    i2_only_v = find(~used2);

    out.idx1_only = orig1(i1_only_v);
    out.idx2_only = orig2(i2_only_v);

    out.counts = struct( ...
        'n12',     numel(out.idx12_1), ...
        'n1_only', numel(out.idx1_only), ...
        'n2_only', numel(out.idx2_only) );

    if nargout <= 2
        varargout = {idx1, idx2};
    elseif nargout == 3
        varargout = {idx1, idx2, d12};
    else
        varargout = {idx1, idx2, d12, out};
    end

    if numel(varargout) < nargout
        varargout(end+1:nargout) = {[]};
    end
    return;
end

% =========================
% 3-channel mode
% =========================
% Strategy:
%   1) one-to-one matches: (1<->2) and (1<->3) via mutual NN
%   2) keep only coord1 points that appear in both match sets
%   3) enforce coord2<->coord3 distance <= r (strict tri-coloc)

[i1_12, i2_12, d12] = mutualNN_withinRadius(c1, c2, r2);
[i1_13, i3_13, d13] = mutualNN_withinRadius(c1, c3, r2);

n1v = size(c1, 1);

% Fast lookup tables keyed by coord1 valid-index
map12_i2  = zeros(n1v, 1);
map12_d12 = NaN(n1v, 1);
map13_i3  = zeros(n1v, 1);
map13_d13 = NaN(n1v, 1);

if ~isempty(i1_12)
    map12_i2(i1_12)  = i2_12;
    map12_d12(i1_12) = d12;
end
if ~isempty(i1_13)
    map13_i3(i1_13)  = i3_13;
    map13_d13(i1_13) = d13;
end

common = (map12_i2 > 0) & (map13_i3 > 0);
i1_common = find(common);

if isempty(i1_common)
    i1_tri = zeros(0,1);
    i2_tri = zeros(0,1);
    i3_tri = zeros(0,1);
    dist123 = zeros(0,3);

    idx1 = zeros(0,1);
    idx2 = zeros(0,1);
    idx3 = zeros(0,1);
else
    i2 = map12_i2(i1_common);
    i3 = map13_i3(i1_common);

    d12c = map12_d12(i1_common);
    d13c = map13_d13(i1_common);

    % Enforce coord2<->coord3 within radius too (strict tri-coloc)
    dx23 = c2(i2,1) - c3(i3,1);
    dy23 = c2(i2,2) - c3(i3,2);
    d23sq = dx23.*dx23 + dy23.*dy23;

    keep = (d23sq <= r2);

    i1_tri = i1_common(keep);
    i2_tri = i2(keep);
    i3_tri = i3(keep);

    d12t = d12c(keep);
    d13t = d13c(keep);
    d23t = sqrt(d23sq(keep));

    idx1 = orig1(i1_tri);
    idx2 = orig2(i2_tri);
    idx3 = orig3(i3_tri);

    dist123 = [d12t, d13t, d23t];
end

% ----- Legacy outputs (ReturnOut=false OR out not requested) -----
if ~wantOut
    if nargout <= 3
        varargout = {idx1, idx2, idx3};
    else
        varargout = {idx1, idx2, idx3, dist123};
    end
    if numel(varargout) < nargout
        varargout(end+1:nargout) = {[]};
    end
    return;
end

% ----- Extended output struct (ReturnOut=true AND requested) -----
% Compute 2<->3 mutual matches (needed for 23-only + singles + exclusions)
[i2_23, i3_23, d23] = mutualNN_withinRadius(c2, c3, r2);

n2v = size(c2,1);
n3v = size(c3,1);

% Build pair maps for combo classification
map21_i1 = zeros(n2v,1);   % coord2 idx -> coord1 idx (from 12 matches)
map23_i3 = zeros(n2v,1);   % coord2 idx -> coord3 idx (from 23 matches)

map31_i1 = zeros(n3v,1);   % coord3 idx -> coord1 idx (from 13 matches)
map32_i2 = zeros(n3v,1);   % coord3 idx -> coord2 idx (from 23 matches)

if ~isempty(i1_12)
    map21_i1(i2_12) = i1_12;
end
if ~isempty(i1_13)
    map31_i1(i3_13) = i1_13;
end
if ~isempty(i2_23)
    map23_i3(i2_23) = i3_23;
    map32_i2(i3_23) = i2_23;
end

% Tri masks per channel (indices into cleaned arrays)
tri1 = false(n1v,1); tri1(i1_tri) = true;
tri2 = false(n2v,1); tri2(i2_tri) = true;
tri3 = false(n3v,1); tri3(i3_tri) = true;

used1 = tri1;
used2 = tri2;
used3 = tri3;

% ---- Pure pair-only buckets (mutually exclusive by construction) ----
% 12-only: pair exists AND neither member has ANY mutual match to ch3 AND not in tri
pure12_mask = false(numel(i1_12),1);
if ~isempty(i1_12)
    pure12_mask = ...
        ~tri1(i1_12) & ~tri2(i2_12) & ...
        (map13_i3(i1_12) == 0) & ...
        (map23_i3(i2_12) == 0);
end
i1_12_only = i1_12(pure12_mask);
i2_12_only = i2_12(pure12_mask);
d12_only   = d12(pure12_mask);

used1(i1_12_only) = true;
used2(i2_12_only) = true;

% 13-only
pure13_mask = false(numel(i1_13),1);
if ~isempty(i1_13)
    pure13_mask = ...
        ~tri1(i1_13) & ~tri3(i3_13) & ...
        (map12_i2(i1_13) == 0) & ...
        (map32_i2(i3_13) == 0);
end
i1_13_only = i1_13(pure13_mask);
i3_13_only = i3_13(pure13_mask);
d13_only   = d13(pure13_mask);

used1(i1_13_only) = true;
used3(i3_13_only) = true;

% 23-only
pure23_mask = false(numel(i2_23),1);
if ~isempty(i2_23)
    pure23_mask = ...
        ~tri2(i2_23) & ~tri3(i3_23) & ...
        (map21_i1(i2_23) == 0) & ...
        (map31_i1(i3_23) == 0);
end
i2_23_only = i2_23(pure23_mask);
i3_23_only = i3_23(pure23_mask);
d23_only   = d23(pure23_mask);

used2(i2_23_only) = true;
used3(i3_23_only) = true;

% ---- Singles (only in one channel) ----
only1_mask = (~used1) & (map12_i2 == 0) & (map13_i3 == 0);
only2_mask = (~used2) & (map21_i1 == 0) & (map23_i3 == 0);
only3_mask = (~used3) & (map31_i1 == 0) & (map32_i2 == 0);

i1_only = find(only1_mask);
i2_only = find(only2_mask);
i3_only = find(only3_mask);

used1(i1_only) = true;
used2(i2_only) = true;
used3(i3_only) = true;

% ---- Excluded conflicts (multi-linked but not classifiable without an "ambiguous" bucket) ----
conf1_mask = ~used1;
conf2_mask = ~used2;
conf3_mask = ~used3;

% Pack out (indices are into ORIGINAL coord arrays, i.e., pre-cleaning rows)
out = struct();
out.mode = '3ch';
out.r    = r;

% strict tri (same as legacy)
out.idx123_1 = idx1(:);
out.idx123_2 = idx2(:);
out.idx123_3 = idx3(:);
out.dist123  = dist123;

% 12-only
out.idx12_only_1 = orig1(i1_12_only);
out.idx12_only_2 = orig2(i2_12_only);
out.dist12_only  = d12_only(:);

% 13-only
out.idx13_only_1 = orig1(i1_13_only);
out.idx13_only_3 = orig3(i3_13_only);
out.dist13_only  = d13_only(:);

% 23-only
out.idx23_only_2 = orig2(i2_23_only);
out.idx23_only_3 = orig3(i3_23_only);
out.dist23_only  = d23_only(:);

% singles
out.idx1_only = orig1(i1_only);
out.idx2_only = orig2(i2_only);
out.idx3_only = orig3(i3_only);

% counts
out.counts = struct( ...
    'n123',     numel(out.idx123_1), ...
    'n12_only', numel(out.idx12_only_1), ...
    'n13_only', numel(out.idx13_only_1), ...
    'n23_only', numel(out.idx23_only_2), ...
    'n1_only',  numel(out.idx1_only), ...
    'n2_only',  numel(out.idx2_only), ...
    'n3_only',  numel(out.idx3_only) );

% excluded (counts only; no ambiguous indices returned)
out.excluded = struct( ...
    'n1', nnz(conf1_mask), ...
    'n2', nnz(conf2_mask), ...
    'n3', nnz(conf3_mask) );

% Return outputs
if nargout <= 3
    varargout = {idx1, idx2, idx3};
elseif nargout == 4
    varargout = {idx1, idx2, idx3, dist123};
else
    varargout = {idx1, idx2, idx3, dist123, out};
end

if numel(varargout) < nargout
    varargout(end+1:nargout) = {[]};
end

end


% =========================
% Helpers
% =========================

function [c, orig] = cleanCoords(coord)
if isempty(coord)
    c    = zeros(0,2);
    orig = zeros(0,1);
    return;
end
if ~isnumeric(coord) || size(coord,2) < 2
    error('Each coord input must be numeric with at least 2 columns: [X Y].');
end
c = double(coord(:,1:2));
valid = isfinite(c(:,1)) & isfinite(c(:,2));
orig = find(valid);
c = c(valid,:);
end


function [iA, iB, dist] = mutualNN_withinRadius(A, B, r2)
% Returns one-to-one mutual nearest neighbor pairs within squared radius r2.
% A: [nA x 2], B: [nB x 2]
% Outputs iA, iB indices into A and B (valid-filtered), and dist (pixels).

nA = size(A,1);
nB = size(B,1);

if nA == 0 || nB == 0
    iA = zeros(0,1);
    iB = zeros(0,1);
    dist = zeros(0,1);
    return;
end

% Fast path: KD-tree nearest neighbor if available
if exist('knnsearch', 'file') == 2
    % Nearest B for each A
    [j_star, dA] = knnsearch(B, A, 'K', 1);
    % Nearest A for each B
    [i_star, ~]  = knnsearch(A, B, 'K', 1);

    rows = (1:nA).';
    mutual = (i_star(j_star) == rows);
    keep = mutual & ((dA .* dA) <= r2);

    iA = rows(keep);
    iB = j_star(keep);
    dist = dA(keep);
    return;
end

% Fallback: blockwise brute force (no toolbox dependency)
X2 = B(:,1).';
Y2 = B(:,2).';

minD2_row = inf(nA,1);
j_star    = ones(nA,1);

minD2_col = inf(1,nB);
i_star    = ones(1,nB);

% Block size to avoid blowing memory on huge nA*nB
maxBlockElems = 1e7;
blockSize = max(1, floor(maxBlockElems / nB));

for s = 1:blockSize:nA
    e = min(nA, s + blockSize - 1);

    X1 = A(s:e,1);
    Y1 = A(s:e,2);

    dx = bsxfun(@minus, X1, X2);
    dy = bsxfun(@minus, Y1, Y2);

    D2 = dx.*dx + dy.*dy;

    % Row minima (A -> B)
    [rowMin, rowArg] = min(D2, [], 2);
    minD2_row(s:e) = rowMin;
    j_star(s:e)    = rowArg;

    % Column minima (B -> A)
    [colMin, colArg] = min(D2, [], 1);
    upd = colMin < minD2_col;
    if any(upd)
        minD2_col(upd) = colMin(upd);
        i_star(upd)    = (s - 1) + colArg(upd);
    end
end

rows = (1:nA).';
mutual = (i_star(j_star).') == rows;
keep = mutual & (minD2_row <= r2);

iA = rows(keep);
iB = j_star(keep);
dist = sqrt(minD2_row(keep));
end
