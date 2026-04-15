function [selected_tracksCoordinates, track_len] = traj_selection_v03_by_length(tracksCoordinates, traj_longer_or_shorter, traj_len_criterion)
% traj_selection_v02_by_length
% Select trajectories by length.
%
% INPUTS
%   tracksCoordinates.X : [nFrames x nTracks] (NaN allowed)
%   tracksCoordinates.Y : [nFrames x nTracks]
%   traj_longer_or_shorter : 1 => keep len >= criterion, 0 => keep len < criterion
%   traj_len_criterion : integer-ish length threshold (frames)
%
% OUTPUTS
%   selected_tracksCoordinates : struct with .X and .Y subset
%   track_len : [nTracks x 1] length per track (#frames with valid X and Y)

% ---- input checks ----
if ~isstruct(tracksCoordinates) || ~isfield(tracksCoordinates,'X') || ~isfield(tracksCoordinates,'Y')
    error('tracksCoordinates must be a struct with fields X and Y.');
end
X = tracksCoordinates.X;
Y = tracksCoordinates.Y;
if ~ismatrix(X) || ~ismatrix(Y) || ~isequal(size(X), size(Y))
    error('tracksCoordinates.X and .Y must be 2D matrices of the same size.');
end

if ~(traj_longer_or_shorter==0 || traj_longer_or_shorter==1)
    error('traj_longer_or_shorter must be 0 or 1.');
end

if ~isscalar(traj_len_criterion) || ~isfinite(traj_len_criterion)
    error('traj_len_criterion must be a finite scalar.');
end
traj_len_criterion = round(traj_len_criterion);

% ---- compute track lengths (vectorized, fast) ----
valid = isfinite(X) & isfinite(Y);
track_len = sum(valid, 1).';   % [nTracks x 1]

% ---- select indices ----
if traj_longer_or_shorter == 1
    idx = (track_len >= traj_len_criterion);
else
    idx = (track_len < traj_len_criterion);
end

% ---- output subset ----
selected_tracksCoordinates = struct();
selected_tracksCoordinates.X = X(:, idx);
selected_tracksCoordinates.Y = Y(:, idx);

% If you have other fields in tracksCoordinates you want to carry along, do it here.
% Example:
% f = fieldnames(tracksCoordinates);
% for k = 1:numel(f)
%     name = f{k};
%     if ~ismember(name, {'X','Y'})
%         selected_tracksCoordinates.(name) = tracksCoordinates.(name);
%     end
% end

end