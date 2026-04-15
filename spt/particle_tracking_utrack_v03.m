function [selected_tracksCoordinates,track_features,processingInfo] = particle_tracking_utrack_v03(processingInfo)
%% Metadata parsing
metadata.voxel_size = processingInfo.AppInputSettings.VoxelSize;
metadata.dt = processingInfo.AppInputSettings.FrameInterval;
metadata.min_track_len = 8;
metadata.filepath = processingInfo.sourceImagePath;
metadata.filenameBase = processingInfo.imageName;
metadata.lastImageNum = processingInfo.nFramesUsed;
metadata.imgWidth = processingInfo.sizeX;
metadata.imgHeight = processingInfo.sizeY;

%% MovieInfo parsing
movieInfo = coords_to_movieInfo_utrack(processingInfo.Coords);

%% Single particle tracking
tic
disp(['Start particle tracking for ' metadata.filenameBase])

% Divide the track to smaller tracks
% Divide movieInfo into multiple cycles --> movieInfo_cycle
track_length = metadata.lastImageNum;   % Default value = metadata.lastImageNum
subtrack_length = 200;                  % Default value = 200 frames
number_of_subtrack = floor(track_length/subtrack_length);
lastsubtrack_firstframe = (number_of_subtrack*subtrack_length)+1;
lastsubtrack_lastframe = track_length;
lastsubtrack_length = lastsubtrack_lastframe - lastsubtrack_firstframe;

if track_length < subtrack_length
    firstframe = 1;
    lastframe = track_length;
    framelist = [firstframe lastframe];
elseif lastsubtrack_length < 5  % If remainder <5 frames --> do not include them
    framelist = [];
    for i = 1:number_of_subtrack
        firstframe = ((i-1)*subtrack_length)+1;
        lastframe = i*subtrack_length;
        A = [firstframe lastframe];
        framelist = [framelist;A];
    end
else    % If remainder >=5 frames --> put them as last track
    framelist = [];
    for i = 1:number_of_subtrack
        firstframe = ((i-1)*subtrack_length)+1;
        lastframe = i*subtrack_length;
        A = [firstframe lastframe];
        framelist = [framelist;A];
    end
    lastsubtrack = [lastsubtrack_firstframe lastsubtrack_lastframe];
    framelist = [framelist;lastsubtrack];
end

% ---- ADD: 10-frame junction windows across subtrack boundaries ----
% For each boundary between consecutive rows in framelist, add a 10-frame window that
% straddles the cut: [boundary-4 ... boundary] U [boundary+1 ... boundary+5]
junction_window = 10;
B_include = ceil(junction_window/2);      % frames on the left (including boundary) -> 5
A_right   = junction_window - B_include;  % frames strictly on the right          -> 5

junctionlist = [];
if size(framelist,1) > 1
    junctionlist = zeros(size(framelist,1)-1, 2);
    for k = 1:size(framelist,1)-1
        boundary = framelist(k,2);              % last frame of subtrack k
        s = boundary - (B_include - 1);         % include the boundary on the left side
        e = boundary + A_right;                 % include 5 frames on the right side

        % Clamp to [1, track_length]
        s = max(1, s);
        e = min(track_length, e);

        % Enforce exact length (10) when possible after clamping
        curLen = e - s + 1;
        if curLen < junction_window
            deficit = junction_window - curLen;
            % Prefer expanding left first (keeps boundary near center), then right
            shiftLeft = min(s - 1, deficit);
            s = s - shiftLeft;
            curLen = e - s + 1;
            if curLen < junction_window
                e = min(track_length, e + (junction_window - curLen));
            end
        end

        junctionlist(k,:) = [s e];
    end

    % Indices of appended junction windows (useful for later linking)
    junction_rows = (size(framelist,1)+1) : (size(framelist,1)+size(junctionlist,1));

    % Append junction windows to framelist
    framelist = [framelist; junctionlist];
else
    junction_rows = [];  % no junctions when only one subtrack
end

% Tracking
parfor i = 1:length(framelist(:,1))
    tracksFinal_divide{i,1} = scriptTrackGeneral_func(movieInfo(framelist(i,1):framelist(i,2)));
end

% Join divided tracks
search_radius = 2;      % Default value = 5 pixels
tracksCoordinates = tracksFinal_conquer_v04(tracksFinal_divide, framelist, junction_rows, search_radius);

%% Calculate track features and remove bg track
track_features = track_feature_v01(tracksCoordinates, metadata);

%% Remove bg track
% by minorAxis 
idx = find(log10(track_features.minorAxis)>-5);
tracksCoordinates.X = tracksCoordinates.X(:,idx);
tracksCoordinates.Y = tracksCoordinates.Y(:,idx);

% by track len
if metadata.min_track_len > 1
    selected_tracksCoordinates = traj_selection_v03_by_length(tracksCoordinates,1,metadata.min_track_len);

    if isempty(selected_tracksCoordinates.X)
        disp(['There is no track longer than ' num2str(metadata.min_track_len) ' frames']);
        selected_tracksCoordinates = tracksCoordinates;
        track_features = track_feature_v01(tracksCoordinates, metadata);
    else
        track_features = track_feature_v01(selected_tracksCoordinates, metadata);
    end
end

%% Save data
processingInfo.min_track_len = metadata.min_track_len;
processingInfo.tracksCoordinates = selected_tracksCoordinates;
processingInfo.track_features = track_features;

elapsedTime = toc;
disp(['Particle tracking for ' metadata.filenameBase ' is completed in ' num2str(elapsedTime,'%.2f') 's'])
