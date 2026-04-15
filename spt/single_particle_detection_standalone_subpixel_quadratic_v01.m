function locs = single_particle_detection_standalone_subpixel_quadratic_v01(movie, centers, varargin)
%SINGLE_PARTICLE_DETECTION_SUBPIXEL_QUADRATIC_V01
% Fast subpixel estimator using quadratic peak interpolation (3x3). This is
% not as precise as single_particle_detection_standalone_subpixel_precision
% code but it run 100 times faster.
%
% locs = single_particle_detection_subpixel_quadratic_v01(movie, centers, ...)
%
% Inputs
%   movie   : HxW single frame OR HxWxT stack
%   centers : Nx3 or Nx4 numeric array: [frame, x_int, y_int, (optional) net_gradient]
%             MATLAB 1-indexed, x=column, y=row
%
% Name-Value options
%   'Clamp'     : clamp subpixel offsets to [-Clamp, +Clamp]. default 0.5
%   'BgMethod'  : 'none' | 'min' | 'medianCorners'. default 'medianCorners'
%   'Use2D'     : true/false, use 2D quadratic fit. default true
%
% Output (table)
%   locs.frame (uint32)
%   locs.x, locs.y (single)      subpixel coordinates (global, 1-indexed)
%   locs.dx, locs.dy (single)    subpixel offsets relative to integer center
%   locs.net_gradient (single)   passed through if provided, else NaN
%   locs.success (logical)       true if a valid estimate was produced

p = inputParser;
p.addParameter('Clamp', 0.5);
p.addParameter('BgMethod', 'medianCorners');
p.addParameter('Use2D', true);
p.parse(varargin{:});
opt = p.Results;

clampVal = single(opt.Clamp);
bgMethod = lower(string(opt.BgMethod));
use2D = logical(opt.Use2D);

if ndims(movie) == 2
    movie = reshape(movie, size(movie,1), size(movie,2), 1);
end
movie = single(movie);
[H,W,T] = size(movie);

if ~isnumeric(centers) || size(centers,2) < 3
    error('centers must be Nx3 or Nx4: [frame x y (net_gradient)].');
end

N = size(centers,1);
frame = int32(centers(:,1));
xInt  = int32(centers(:,2));
yInt  = int32(centers(:,3));

if size(centers,2) >= 4
    netg = single(centers(:,4));
else
    netg = single(nan(N,1));
end

% Valid indices for 3x3 neighborhood
valid = frame >= 1 & frame <= T & xInt >= 2 & xInt <= (W-1) & yInt >= 2 & yInt <= (H-1);

% Prepare outputs
xSub = single(nan(N,1));
ySub = single(nan(N,1));
dx   = single(nan(N,1));
dy   = single(nan(N,1));
success = false(N,1);

if ~any(valid)
    locs = table(uint32(frame), xSub, ySub, dx, dy, netg, success, ...
        'VariableNames', {'frame','x','y','dx','dy','net_gradient','success'});
    return;
end

% Flatten movie for fast gather
M = movie(:);
HW = int32(H*W);

f = frame(valid);
x = xInt(valid);
y = yInt(valid);

% Linear index (MATLAB column-major):
% idx = y + (x-1)*H + (t-1)*H*W
idx = y + (x-1)*int32(H) + (f-1)*HW;

% Offsets for 3x3 patch in this order:
% 1 TL, 2 T, 3 TR, 4 L, 5 C, 6 R, 7 BL, 8 B, 9 BR
off = int32([-1-int32(H), -1, -1+int32(H), ...
    -int32(H), 0, int32(H), ...
    1-int32(H), 1, 1+int32(H)]);

Z = single(zeros(numel(idx), 9));
for j = 1:9
    Z(:,j) = M(idx + off(j));
end

% Optional background subtraction (very cheap)
switch bgMethod
    case "none"
        bg = single(0);
    case "min"
        bg = min(Z, [], 2);
        Z = Z - bg;
    case "mediancorners"
        % corners: TL, TR, BL, BR => cols 1,3,7,9
        corners = sort(Z(:,[1 3 7 9]), 2);
        bg = 0.5*(corners(:,2) + corners(:,3)); % median of 4
        Z = Z - bg;
    otherwise
        error("BgMethod must be 'none', 'min', or 'medianCorners'.");
end

% Keep nonnegative (helps stability when bg subtraction overshoots)
Z(Z < 0) = 0;

% 1D quadratic offsets (fallback / also useful as sanity)
I0 = Z(:,5);
IL = Z(:,4); IR = Z(:,6);
IT = Z(:,2); IB = Z(:,8);

denx = (IL - 2*I0 + IR);
deny = (IT - 2*I0 + IB);

dx1 = (IL - IR) ./ (2*denx);
dy1 = (IT - IB) ./ (2*deny);

% Valid maxima curvature (den < 0) and finite
ok1x = isfinite(dx1) & denx < 0;
ok1y = isfinite(dy1) & deny < 0;
dx1(~ok1x) = 0;
dy1(~ok1y) = 0;

% Clamp
dx1 = max(min(dx1, clampVal), -clampVal);
dy1 = max(min(dy1, clampVal), -clampVal);

dxv = dx1;
dyv = dy1;
ok = ok1x & ok1y;

if use2D
    % 2D quadratic fit z = ax^2 + by^2 + cxy + dx + ey + f
    % Precomputed least-squares projection for the fixed 3x3 grid.
    xv = [-1 0 1 -1 0 1 -1 0 1]';
    yv = [-1 -1 -1 0 0 0 1 1 1]';
    A = [xv.^2, yv.^2, xv.*yv, xv, yv, ones(9,1)];
    B = (A.'*A)\A.';            % 6x9
    coef = double(Z) * B.';     % Nx6

    a = coef(:,1); b = coef(:,2); c = coef(:,3);
    d = coef(:,4); e = coef(:,5);

    A2 = 2*a; B2 = 2*b; C = c;
    detH = A2.*B2 - C.^2;

    % Vertex solve:
    % [A2 C; C B2] [x0;y0] = [-d; -e]
    x0 = (-d.*B2 + C.*e) ./ detH;
    y0 = (-A2.*e + C.*d) ./ detH;

    % For a maximum, Hessian should be negative definite: A2<0, B2<0, det>0
    ok2 = isfinite(x0) & isfinite(y0) & (A2 < 0) & (B2 < 0) & (detH > 1e-12);

    % Clamp offsets (keeps nonsense from poisoning you)
    x0 = max(min(single(x0), clampVal), -clampVal);
    y0 = max(min(single(y0), clampVal), -clampVal);

    % Use 2D where valid, else fallback to 1D
    dxv(ok2) = x0(ok2);
    dyv(ok2) = y0(ok2);
    ok = ok2 | ok; % success if either method produced something sane
end

% Write back into full-size outputs
idxValid = find(valid);
dx(idxValid) = dxv;
dy(idxValid) = dyv;
xSub(idxValid) = single(x) + dxv;
ySub(idxValid) = single(y) + dyv;
success(idxValid) = ok;

locs = table(uint32(frame), xSub, ySub, dx, dy, netg, success, ...
    'VariableNames', {'frame','x','y','dx','dy','net_gradient','success'});
end