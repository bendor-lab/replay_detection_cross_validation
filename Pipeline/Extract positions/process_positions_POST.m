function process_positions_POST(track_pairs)
% Inputs:
% loads 'position' (structure created in process_positions_PRE.m), 'positions' (structure with raw x,y,t extracted from video)
% and list_of_parameters.m
%track_pairs : 'If there is re-exposure to the same track, please write which tracks are the same (e.g. [1,3; 2,4]). Else, input [].'
% Output: 'position' structure containing new x,y,v and linearized position

parameters= list_of_parameters;
max_distance= parameters.max_pixel_distance;  % that is max distance from track mask
max_distance_jump= parameters.max_pixel_distance_jump; % for brief jumps away
position_filter_length= parameters.position_filter_length; % to smooth speed data

load('position_data');
load('extracted_video_data');

number_of_tracks= length(position.linear);
[lin_timestamps_session, unique_idx]= unique(position.clean.t); %unique_indx tells you if there's repeated timestamps
if length(lin_timestamps_session)~=length(position.clean.t)
    disp('ERROR- repeated timestamps')
end

for r=1:number_of_tracks
    % Clean variable 'linear' to get rid of previous track
    linear=[];
    
    % find indices on track
    track_indices= position.linear(r).cropped_indices;
    x_bins=1:ceil(max(position.clean.x(track_indices)));
    y_bins=1:ceil(max(position.clean.y(track_indices)));
    
    % Use cleaned track data and find closest point to template track for each
    % case. Get output structure called grid
    T= position.clean.t(track_indices);
    X= position.clean.x(track_indices);
    Y= position.clean.y(track_indices);
    grid(r).index=zeros(length(x_bins),length(y_bins));
    grid(r).distance=zeros(length(x_bins),length(y_bins));
    
    % If there's been re-exposure, the mask from second exposure will be used for the first linearizing the first exposure.
    if ~isempty(track_pairs)
        [track_indx,~] = find(r == track_pairs); % find the row indx
        track_mask_x = position.linear(track_pairs(track_indx,2)).track_mask_x;
        track_mask_y = position.linear(track_pairs(track_indx,2)).track_mask_y;
    else
        track_mask_x = position.linear(r).track_mask_x;
        track_mask_y = position.linear(r).track_mask_y;
    end
    %calculate distance from track
    for i=1:length(x_bins)
        for j=1:length(y_bins)
            
            distance = sqrt((x_bins(i)-track_mask_x).^2+(y_bins(j)-track_mask_y).^2);
            distance_score = max_distance-distance;
            distance_score((distance_score<0))=0; %distance score is maximum when 0 distance from maze
            if isempty(distance_score)
                grid(r).index(i,j) = NaN;
                grid(r).distance(i,j) = NaN;
            else  % take weighted average (index*distance_score) to avoid jumps around corners of maze
                grid(r).index(i,j) = sum((1:length(distance_score)).*distance_score)/sum(distance_score);
                grid(r).distance(i,j) = mean(distance_score((distance_score>0)));
            end
        end
    end
    position.linear(r).grid= grid(r);
    
    % Matrix with equivalence between interpolated timestamps,x,y and new pixels(from linear.grid)
    for j=1:length(T)   %finds the equivalence between a pixel and 'x and y' position
        linear.x(j)=grid(r).index(round(X(j)),round(Y(j)));
        linear.d(j)=grid(r).distance(round(X(j)),round(Y(j)));
    end
    
    % remove NaNs (not on track)
    index_nan= find(isnan(linear.x)==1);
    if ~isempty(index_nan)
        linear.x(index_nan)=[];
        linear.d(index_nan)=[];
        track_indices(index_nan)=[];
        T(index_nan)=[];
        X(index_nan)=[];
        Y(index_nan)=[];
    end
    
    % eliminate large jumps: don't put too small a value for max_distance_jump
    raw_lin_timestamps{r}= T;   %original values before removing jumps
    raw_lin_positions{r}= linear.x;
    lin_timestamps{r}= T;
    lin_positions{r}= linear.x;
    x_positions{r}= X;
    y_positions{r}= Y;
    number_of_jumps_removed= 0;  % initialise
    jumps_removed_idx{r}=[];
    jumps= find(abs(diff(lin_positions{r})) > max_distance_jump);
    jumps_to_remove_index=[];
    
    for ii = 1:length(jumps)
        this_jump = jumps(ii);
        if this_jump == 1 %if the jump is between the first two points
            jumps_to_remove_index=[jumps_to_remove_index this_jump];
        elseif this_jump >= (length(lin_positions{r})-1)
            jumps_to_remove_index=[jumps_to_remove_index this_jump];
        elseif abs(lin_positions{r}(this_jump)- lin_positions{r}(this_jump+2))< max_distance_jump %if it's a one point jump
            jumps_to_remove_index=[jumps_to_remove_index this_jump+1];
            jumps_removed_idx{r}= [jumps_removed_idx{r}, this_jump+1];
        end
    end
    
    lin_timestamps{r}(jumps_to_remove_index)=[];
    lin_positions{r}(jumps_to_remove_index) =[];
    x_positions{r}(jumps_to_remove_index)   =[];
    y_positions{r}(jumps_to_remove_index)   =[];
    track_indices(jumps_to_remove_index)    =[];
    
    % Second method to remove jumps
    c=1; jumps=[];
    for i=7:(length(lin_timestamps{r})-6)
        if abs(median(lin_positions{r}((i-6):(i+6)))-lin_positions{r}(i))> parameters.max_distance_jump
            jumps(c)=i;
            c=c+1;
        end
    end
    
    % number_of_jumps_removed=length(jumps);
    lin_timestamps{r}(jumps)=[];
    lin_positions{r}(jumps)=[];
    x_positions{r}(jumps)=[];
    y_positions{r}(jumps)=[];
    track_indices(jumps)=[];
    jumps_removed_idx{r}= [jumps_removed_idx{r}, jumps];
    
    % convert linear position to cm by normalize linear position and multiply by length (in cm) 
   lin_positions{r}= 100*position.linear(r).length*(lin_positions{r}-min(lin_positions{r}))/(max(lin_positions{r})-min(lin_positions{r}));

    disp(strcat('Track ',num2str(r), ' number of jumps removed'))
    disp(length(jumps_removed_idx{r}))
end

% Create super regular timestamps using raw timestamps
lin_timestamps_step = 0.04; % our camera is at 25Hz
position.t = min(positions.timestamps):lin_timestamps_step:max(positions.timestamps);%timestamps

% Initialize X&Y positions with length of raw timestamps
X_positions = NaN(size(lin_timestamps_session));
Y_positions = NaN(size(lin_timestamps_session));

f= figure('Name','Interpolated dropped samples vs raw');
% interpolate track_indices to upsampled version
for kk=1:number_of_tracks
    [time_diff,~]= min(abs(lin_timestamps_session-min(lin_timestamps{kk}))); % should be perfectly aligned
    if time_diff ~= 0
        disp('error in synchronisation');
        keyboard;
    end
    [~,~,idx_track]= intersect(lin_timestamps{kk},lin_timestamps_session);
    if length(idx_track)~= length(lin_timestamps{kk})
        disp('Something might have gone wrong during deletion of repeated timestamps. You might want to double-check it')
    end
    
    % Fill X&Y track positions
    X_positions(idx_track) = x_positions{kk};
    Y_positions(idx_track) = y_positions{kk};
    position.linear(kk).linear =interp1(lin_timestamps{kk}, lin_positions{kk},position.t,'nearest'); % interpolate linear track
    
    % Visualise changes in tracking caused by interpolation and removing jumps
    set(0, 'currentfigure', f);
    ax(2*kk-1) = subplot(number_of_tracks,2,2*kk-1);
    plot(position.t,position.linear(kk).linear,'k.');
    yl = ylim;
    title(['Track' num2str(kk) ' interpolated']);
    
    ax(2*kk) = subplot(number_of_tracks,2,2*kk);
    raw_lin_positions{kk}=  position.linear(r).cm_per_pixel*(raw_lin_positions{kk}-min(lin_positions{kk}));
    plot(raw_lin_timestamps{kk},raw_lin_positions{kk},'k.');
    hold on;
    plot(raw_lin_timestamps{kk}(jumps_removed_idx{kk}),raw_lin_positions{kk}(jumps_removed_idx{kk}),'rx');
    ylim([yl(1) yl(2)+20]);
    legend('raw','dropped');
    title(['Max jump: '  num2str(max_distance_jump) ', Percent dropped: ' num2str(100*length(jumps_removed_idx{kk})/length(raw_lin_timestamps{kk})) '%']);
    linkaxes(ax(2*kk-1:2*kk),'x');
end

% Interpolate sleepbox indices to upsampled version & align to main timescale
sleepbox= position.sleep.t;
[time_diff,~]= min(abs(lin_timestamps_session-min(sleepbox))); % should be perfectly aligned
if time_diff ~= 0
    disp('error in synchronisation');
    keyboard;
end
[~,~,idx_sleep]= intersect(sleepbox,lin_timestamps_session);
sleepbox_temp= NaN(size(lin_timestamps_session));
sleepbox_temp(idx_sleep)= sleepbox;
position.sleepbox= interp1(lin_timestamps_session,sleepbox_temp,position.t,'nearest');
X_positions(idx_sleep) = position.clean.x(idx_sleep);
Y_positions(idx_sleep) = position.clean.y(idx_sleep);

% Remove NaNs from non desirable points
XY_nan= (isnan(X_positions) | isnan(Y_positions));
X_positions(XY_nan)=[];
Y_positions(XY_nan)=[];
XY_timestamps_session= lin_timestamps_session;
XY_timestamps_session(XY_nan)=[];

% save new track indices and track timestamps
for track_id=1:number_of_tracks
    position.linear(track_id).timestamps = position.t(~isnan(position.linear(track_id).linear));
    position.linear(track_id).clean_track_Indices = find(~isnan(position.linear(track_id).linear));
    position.linear(track_id).clean_track_LogicalIndices = ~isnan(position.linear(track_id).linear);
end



% Correct tracking noise inside sleep box
vraw  = sqrt(diff(X_positions').^2+diff(Y_positions').^2)'./lin_timestamps_step; %calculated from clean X & Y
v_cmraw = vraw*mean([position.linear(:).cm_per_pixel]); % convert to cm/s
noise_indices = find(v_cmraw > 300);  % find indices of velocity over threshold
t2 = XY_timestamps_session;
t2(noise_indices) = [];
y2 = Y_positions;
y2(noise_indices) = [];
x2 = X_positions;
x2(noise_indices) = [];

%%% Interpolate and get final clean x, y and velocity
position.x = interp1(t2,x2,position.t,'nearest');
position.y = interp1(t2,y2,position.t,'nearest');
nan_index = isnan(position.x);
if any(isnan(position.x ) | isnan(position.y))
    disp('ERROR- interpolation has produced NaNs- expanding time limits')
    position.x = interp1(t2,x2,position.t,'nearest','extrap'); %extrapolate values outside of range
    position.y = interp1(t2,y2,position.t,'nearest','extrap');
    figure
    plot(position.t,position.x,'b')
    hold on
    plot(position.t,position.y,'r')
    plot(position.t(nan_index),position.x(nan_index),'gx')
    plot(position.t(nan_index),position.y(nan_index),'gx')
    title('x indicates where interpolation produced NaNs before expanding time limits')
end

position.v = sqrt(diff(smooth(position.x',parameters.position_filter_length)).^2+diff(smooth(position.y',parameters.position_filter_length)).^2)'./lin_timestamps_step; %calculated from smoothed clean X & Y
position.v = interp1([position.t(1)-lin_timestamps_step/2 position.t+lin_timestamps_step/2],[position.v(1) position.v position.v(end)],position.t,'nearest');
position.v_cm = position.v*mean([position.linear(:).cm_per_pixel]);
position.cm_per_pixel = mean([position.linear(:).cm_per_pixel]);
position.v_zscore = zscore(position.v);
position.jumps_removed_idx = jumps_removed_idx;

save('extracted_position','position');
end