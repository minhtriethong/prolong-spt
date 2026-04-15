function centers = single_particle_detection_standalone_v01(movie, boxSize, minNetGradient, minIntensity)
% single particle detection using local maxima + optional net-gradient filter.
%
% centers = single_particle_detection_standalone_v01(movie, boxSize, minNetGradient, minIntensity)
%
% Inputs
%   movie          : HxW single frame OR HxWxT stack
%   boxSize        : odd integer (e.g., 5). Defines the local-max neighborhood.
%   minNetGradient : scalar. If finite, filters candidates using Picasso-style net gradient.
%                    Use -Inf (recommended for max speed) to skip net-gradient computation.
%   minIntensity   : scalar. Optional cheap pre-filter (recommended). Use -Inf to disable.
%
% Output
%   centers : Nx3 array [frame, x, y, netgradient] (1-indexed). x = column, y = row.
%
% Notes
%   - For maximum speed: minNetGradient = -Inf and set a sensible minIntensity.
%   - This intentionally does NOT extract 5x5 patches; it only returns center positions.

if nargin < 2 || isempty(boxSize),        boxSize = 5;          end
if nargin < 3 || isempty(minNetGradient), minNetGradient = -Inf; end
if nargin < 4 || isempty(minIntensity),   minIntensity = -Inf;   end

if mod(boxSize, 2) == 0
    error('boxSize must be an odd integer.');
end
if ndims(movie) == 2
    movie = reshape(movie, size(movie,1), size(movie,2), 1);
end

[H, W, T] = size(movie);
r = floor(boxSize/2);

% Precompute unit vectors (Picasso-style) and linear offsets for patch gathering
[dX, dY] = meshgrid(-r:r, -r:r);          % dX = col offsets, dY = row offsets
uX = -single(dX);                          % vector pointing toward center
uY = -single(dY);
nrm = sqrt(uX.^2 + uY.^2);
nrm(r+1, r+1) = 1;                         % avoid division by zero at center
uX = uX ./ nrm;  uY = uY ./ nrm;
uX(r+1, r+1) = 0; uY(r+1, r+1) = 0;        % center excluded

wX = uX(:);                                % 25x1
wY = uY(:);                                % 25x1
offsets = single(dY + dX*H);               % linear offsets in MATLAB (col-major)
offsets = offsets(:)';                     % 1x(boxSize^2)

centersCell = cell(T, 1);

% Border needed:
% - local max needs r pixels margin for full box
% - net gradient needs an extra 1 pixel margin for central-diff gradients
border = r + double(isfinite(minNetGradient));

for t = 1:T
    I = single(movie(:,:,t));

    % Local maxima in boxSize x boxSize using separable moving max (fast, base MATLAB)
    localMax = movmax(movmax(I, [r r], 1), [r r], 2);

    cand = (I == localMax) & (I >= minIntensity);

    % Remove borders so the neighborhood and (optional) gradients are valid
    if border > 0
        cand(1:border, :) = false;
        cand(end-border+1:end, :) = false;
        cand(:, 1:border) = false;
        cand(:, end-border+1:end) = false;
    end

    if ~any(cand(:))
        centersCell{t} = zeros(0,3);
        continue;
    end

    [y, x] = find(cand); % y=row, x=col

    % Central-difference gradients (Picasso-style)
    gY = zeros(H, W, 'single');
    gX = zeros(H, W, 'single');
    gY(2:end-1, :)   = I(3:end, :)   - I(1:end-2, :);
    gX(:, 2:end-1)   = I(:, 3:end)   - I(:, 1:end-2);

    idx = sub2ind([H W], y, x);       % Nx1
    patchIdx = idx + offsets;         % Nx(boxSize^2) via implicit expansion

    % Net gradient per candidate (vectorized)
    ng = gY(patchIdx) * wY + gX(patchIdx) * wX;

    keep = (ng > minNetGradient);
    x = x(keep);
    y = y(keep);

    n = numel(x);
    centersCell{t} = [t*ones(n,1), x, y, ng];
end

centers = vertcat(centersCell{:});
end