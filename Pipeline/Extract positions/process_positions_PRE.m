function process_positions_PRE

% Input: data extracted from NLX 'VT1.nvt'
% - x: x position of largest tracked blob
% - y: y position of largest tracked blob
% - timestamps: timestamps (in seconds)

% Output: linear structure containing new X,Y,T and linearized position template
% this function calls 'crop_in_time.m'

load extracted_video_data
parameters= list_of_parameters;

% plot position on top of extracted frame for reference
if exist('extracted_video_frame.mat')==2
    load extracted_video_frame
    figure('Name','Session');     % plot position on top of extracted frame for reference
    imagesc(extracted_video_frame);
    [~, folder_name]= fileparts(pwd);
    title(['Session ' folder_name]);
    hold on
    plot(positions.x,positions.y,'r.','MarkerSize',1);
    axis xy
else
    figure;
    subplot(1,2,1)
    plot(positions.timestamps,positions.y,'r.','MarkerSize',1); 
    hold on
    plot(positions.timestamps,positions.x,'b.','MarkerSize',1); 
    subplot(1,2,2)
    plot(positions.x,positions.y,'r.','MarkerSize',1);  
end

% this determines how many loops the programme does (one per exposure)
reply = input('How many exposures are there:   ' ,'s');
number_of_tracks = str2num(reply);
out.number_of_tracks = number_of_tracks;

% initialise vectors which will contain cleaned positions (track + sleep)
clean.t = NaN(size(positions.timestamps));
clean.x = NaN(size(positions.x));
clean.y = NaN(size(positions.y));

all_cropped_indices = [];
for r = 1:number_of_tracks
    % select in time each maze in temporal order
    [cropped_indices,cropped_logical] = crop_in_time([positions.x' positions.y'],positions.timestamps,r);
    position.linear(r).cropped_indices = cropped_indices;
    position.linear(r).cropped_logical = cropped_logical;
    all_cropped_indices = [all_cropped_indices cropped_indices];
    
    cropped_x = NaN(size(positions.x));
    cropped_y = NaN(size(positions.y));
    cropped_timestamps = NaN(size(positions.timestamps));
    cropped_x(cropped_indices) = positions.x(cropped_indices); % fill NaNs structure with the cropped values from each track
    cropped_y(cropped_indices) = positions.y(cropped_indices);
    cropped_timestamps(cropped_indices) = positions.timestamps(cropped_indices);
    
    % now select in space
    disp('CROPPING: select polygon of included pixels, click on (-x,-y) portion of plot to finish');
    figure;
    histoplot(cropped_x,cropped_y);
    xlim([0 max(positions.x)]);
    ylim([0 max(positions.y)]);
    annotation('textbox', [0.02 0 0 0.1], 'string', 'Skip/End','LineStyle','none');
    
    new_x=[];  new_y=[];
    xl=[]; yl=[];
    xl=xlim; yl=ylim;
    title('select polygon of included pixels');
    current_point=ginput(1); %change for number of points you want
    while current_point(1)>=min(xl) | current_point(2)>=min(yl)
        new_x=[new_x current_point(1)];
        new_y=[new_y current_point(2)];
        hold on
        plot(new_x,new_y,'r')
        current_point=ginput(1);
    end
    hold on; plot(new_x,new_y,'gx');
    
    if length(new_x)>=3 % if you do not make at least a triangle, take all points
        in = inpolygon(cropped_x,cropped_y,new_x,new_y); % select values inside polygon
        indx_selected_track = find(in==1);
        indx_not_selected_track = find(in==0);
        position.linear(r).cropped_indices = intersect(position.linear(r).cropped_indices,indx_selected_track);
        position.linear(r).cropped_logical(indx_not_selected_track)= 0;
        clean.t(indx_selected_track) = positions.timestamps(indx_selected_track);
        clean.x(indx_selected_track) = positions.x(indx_selected_track);
        clean.y(indx_selected_track) = positions.y(indx_selected_track);
    else % if you have skipped creating a polygon
        clean.t(cropped_indices) = positions.timestamps(cropped_indices);
        clean.x(cropped_indices) = positions.x(cropped_indices);
        clean.y(cropped_indices) = positions.y(cropped_indices);
        idx_zero= find(clean.x==0 | clean.y==0);  %remove zeros
        clean.t(idx_zero)= NaN;
        clean.x(idx_zero)= NaN;
        clean.y(idx_zero)= NaN;
    end
    
    % Re-crop in space if necessary
    data_matrix= [clean.x', clean.y'];
    for i=1:size(data_matrix,2)
        figure;
        plot(positions.timestamps(position.linear(r).cropped_indices),data_matrix(position.linear(r).cropped_indices,i));
        hold on;
        max_trace= max(data_matrix(:,i))+ 5;
        ylim([0 max(data_matrix(:,2) + max_trace)])
        annotation('textbox', [0.01 0 0 0.1], 'string', 'Skip/End x2','LineStyle','none');
        new_x4=[];  new_y4=[];
        xl=[]; yl=[];
        xl=xlim; yl=ylim;
        title('select polygon of included pixels');
        current_point=ginput(1); %change for number of points you want
        while current_point(1)>=min(xl) | current_point(2)>=min(yl) %
            new_x4=[new_x4 current_point(1)];
            new_y4=[new_y4 current_point(2)];
            hold on
            plot(new_x4,new_y4,'r')
            current_point=ginput(1);
        end
        hold on; plot(new_x4,new_y4,'gx')
        
        if length(new_x4)>=3 % if you do not make at least a triangle, take all points
            in= inpolygon(positions.timestamps,data_matrix(:,i),new_x4,new_y4);
            selected_track{i}= find(in==1);
            not_selected_track{i}= find(in==0);
        else  %if skipped selecting a polygon
            selected_track{i}= position.linear(r).cropped_indices;
            not_selected_track{i}= [];
        end
    end
    
    indx_selected_track= intersect(selected_track{1},selected_track{2});  %find insection of indicies in x and y position
    indx_not_selected_track = 1:length(positions.timestamps);
    indx_not_selected_track(indx_selected_track)=[];
    position.linear(r).cropped_indices = intersect(position.linear(r).cropped_indices,indx_selected_track);
    position.linear(r).cropped_logical(indx_not_selected_track)= 0;
    
    clean.t(indx_selected_track)= positions.timestamps(indx_selected_track);
    clean.x(indx_selected_track)= positions.x(indx_selected_track);
    clean.y(indx_selected_track)= positions.y(indx_selected_track);
    position.clean= clean;
    figure
    plot(position.clean.x(indx_selected_track),position.clean.y(indx_selected_track),'r.');
    hold on
    plot(position.clean.x(indx_selected_track),position.clean.y(indx_selected_track),'k');

    %%%%%%%%%% Linearize tracks
    figure;
    histoplot(clean.x(indx_selected_track),clean.y(indx_selected_track));
    xlim([0 max(positions.x)]);
    ylim([0 max(positions.y)]);
    title('linearize track');
    annotation('textbox', [0.02 0 0 0.1], 'string', 'Skip/End x2','LineStyle','none');
    
    new_x2=[];  new_y2=[];
    xl=[]; yl=[];
    xl=xlim; yl=ylim;
    current_point=ginput(1); %change for number of points you want
    while current_point(1)>=min(xl)  & current_point(2)>=min(yl)
        new_x2=[new_x2 current_point(1)];
        new_y2=[new_y2 current_point(2)];
        hold on
        plot(new_x2,new_y2,'r')
        current_point=ginput(1);
    end
    if length(new_x2)<2 % not linearised anything...
        disp('error, you have not linearised the maze');
        keyboard;
    end
    % this is if you want to linearise a multi-arm maze
    branch_type = input(strcat(['is it a multi-arm maze? [y/n]   ']) ,'s');
    if strcmp(branch_type,'y')  %if multi-arm maze
        if current_point(1) <=0
            disp( 'this is if you want to linearise a multi-arm maze');
            new_x3=[];  new_y3=[];
            current_point=ginput(1); %change for number of points you want
            while current_point(1)>=0 && current_point(2)>=0
                new_x3=[new_x3 current_point(1)];
                new_y3=[new_y3 current_point(2)];
                hold on
                plot(new_x3,new_y3,'b')
                current_point=ginput(1);
            end
        end
    end
    % Creates new x and y values
    lin_track(r).x = [];
    lin_track(r).y = [];
    pixel_distance = [];
    pixel_distance_seg1 = []; pixel_distance_seg2 = [];
    for i=1:(length(new_x2)-1)
        max_x = ceil(abs(new_x2(i)-new_x2(i+1)));  %round towards plus infinity
        max_y = ceil(abs(new_y2(i)-new_y2(i+1)));
        steps = round(sqrt(max_x.^2+ max_y.^2));
        lin_track(r).x = [lin_track(r).x linspace(round(new_x2(i)),round(new_x2(i+1)),steps)]; %steps are number of points between start and end of vector
        lin_track(r).y = [lin_track(r).y linspace(round(new_y2(i)),round(new_y2(i+1)),steps)];
        pixel_distance_seg1(i) = sqrt((new_x2(i)-new_x2(i+1)).^2+(new_y2(i)-new_y2(i+1)).^2);
    end
    % idem: runs for multi-arm only
    if strcmp(branch_type,'y')  %if multi-arm maze
        for i=1:(length(new_x3)-1)
            max_x = ceil(abs(new_x3(i)-new_x3(i+1)));
            max_y = ceil(abs(new_y3(i)-new_y3(i+1)));
            steps = round(sqrt(max_x.^2+ max_y.^2));
            lin_track(r).x = [lin_track(r).x linspace(round(new_x3(i)),round(new_x3(i+1)),steps)];
            lin_track(r).y = [lin_track(r).y linspace(round(new_y3(i)),round(new_y3(i+1)),steps)];
            pixel_distance_seg2(i) = sqrt((new_x3(i)-new_x3(i+1)).^2+(new_y3(i)-new_y3(i+1)).^2);
        end
    end
    %remove potential superimposed values in multi-arm case
    [~,idx_unique,~]= unique([lin_track(r).x;lin_track(r).y]','rows','stable'); % trying to remove potential superimposed values in multi-arm case
    position.linear(r).track_mask_x= lin_track(r).x(idx_unique);
    position.linear(r).track_mask_y= lin_track(r).y(idx_unique);
    hold on; plot(lin_track(r).x,lin_track(r).y,'gx')
    
    % get maze length and conversion to px/cm
    switch branch_type  %branch-type is 'y' for a multi-arm maze
        case 'n'
            position.linear(r).branch_type= 'mono';
            track_length = input(strcat(['what is the length of track ' num2str(r) ' in meters:   ']) ,'s');
            position.linear(r).length= str2num(track_length);
            position.linear(r).cm_per_pixel= (100*str2num(track_length))/sum(pixel_distance_seg1); %cm per pixel
        case 'y'
            position.linear(r).branch_type= 'multi_arm';
            seg_length = input(strcat(['what is the length first segment ' num2str(r) ' in meters:   ']) ,'s');
            position.linear(r).cm_per_pixel= (100*str2num(seg_length))/sum(pixel_distance_seg1); %cm per pixel
            track_length = input(strcat(['what is the length of track/arm ' num2str(r) ' in meters:   ']) ,'s');
            position.linear(r).length= str2num(track_length);
    end
    % Close unecessary figures
    fig1h = findall(0,'type','figure','Name','Session');  %Keep this open as its the main GUI figure
    figh = findall(0,'type','figure');
    other_figures = setdiff(figh, fig1h);
    delete(other_figures);
end

%%%%define sleep box
reply = input('How many sleep boxes are there:   ' ,'s');
number_of_sleep_boxes = str2num(reply);
out.number_of_sleep_boxes = number_of_sleep_boxes;

% Replaces track periods by NaNs
sleep.x= positions.x;
sleep.y= positions.y;
sleep.t= positions.timestamps;
sleep.x(all_cropped_indices)= NaN(size(all_cropped_indices));
sleep.y(all_cropped_indices)= NaN(size(all_cropped_indices));
sleep.t(all_cropped_indices)= NaN(size(all_cropped_indices));
position.all_cropped_indices= all_cropped_indices;

figure;
histoplot(sleep.x,sleep.y);
xlim([0 max(positions.x)]);
ylim([0 max(positions.y)]);
annotation('textbox', [0.02 0 0 0.1], 'string', 'Skip/End');
title('tracking data');
indx_selected_sleep=[];
for r=1:number_of_sleep_boxes
    new_x=[]; new_y=[];
    xl=[]; yl=[];
    xl=xlim;
    yl=ylim;
    disp(strcat(['select polygon of sleep box ' num2str(r)]));
    current_point=ginput(1); %change for number of points you want
    while current_point(1)>=min(xl) && current_point(2)>=min(yl)
        new_x=[new_x current_point(1)];
        new_y=[new_y current_point(2)];
        hold on
        plot(new_x,new_y,'r')
        current_point=ginput(1);
    end
    
    hold on; plot(new_x,new_y,'mx')
    in=inpolygon(sleep.x,sleep.y,new_x,new_y);
    indx_selected_sleep= [indx_selected_sleep find(in==1)];
end
    indx_not_selected_sleep= 1:length(sleep.x);
    indx_not_selected_sleep(indx_selected_sleep)=[];
    
    % remove points not selected
    sleep.x(indx_not_selected_sleep)= NaN(size(indx_not_selected_sleep));
    sleep.y(indx_not_selected_sleep)= NaN(size(indx_not_selected_sleep));
    sleep.t(indx_not_selected_sleep)= NaN(size(indx_not_selected_sleep));
    sleep.index=indx_selected_sleep;
    
    %concatenate everything to have both sleep and
    % track positions in variable 'clean'
    temp_concat_x= nansum([clean.x; sleep.x]);
    temp_concat_x(temp_concat_x == 0)= NaN;
    temp_concat_y= nansum([clean.y; sleep.y]);
    temp_concat_y(temp_concat_y == 0)= NaN;
    temp_concat_t= nansum([clean.t; sleep.t]);
    temp_concat_t(temp_concat_t == 0)= NaN;
    
    position.clean.x= temp_concat_x;
    position.clean.y= temp_concat_y;
    position.clean.t= temp_concat_t;
    position.sleep.x= sleep.x;
    position.sleep.y= sleep.y;
    position.sleep.t= sleep.t;
    position.sleep.index= sleep.index;
    position.clean.index= sort([all_cropped_indices sleep.index]);
    
    position.clean.x(isnan(position.clean.t))= [];
    position.clean.y(isnan(position.clean.t))= [];
    position.clean.t(isnan(position.clean.t))= [];
    position.sleep.x(isnan(position.sleep.t))= [];
    position.sleep.y(isnan(position.sleep.t))= [];
    position.sleep.t(isnan(position.sleep.t))= [];
    position.sleep.index(isnan(position.sleep.t))= [];
    position.clean.index(isnan(position.clean.t))= [];
    
    [~,interpolated_indices] = setdiff(positions.timestamps,position.clean.t); %finds indices of the timestamps that will be interpolated
    position.clean.x = interp1(position.clean.t, position.clean.x,positions.timestamps,'linear');
    position.clean.y = interp1(position.clean.t, position.clean.y,positions.timestamps,'linear');
    position.clean.t = positions.timestamps;
    position.interpolated_indices = interpolated_indices; 

close all;
%plot finalized position data that will be saved
figure
number_of_tracks= length(position.linear);
subplot(1,3,1)
hold on
title('x vs y')
plot(position.clean.x,position.clean.y,'k.')
subplot(1,3,2)
hold on
title('t vs x')
plot(position.clean.t,position.clean.x,'k.')
subplot(1,3,3)
hold on
title('t vs y')
plot(position.clean.t,position.clean.y,'k.')
for r=1:number_of_tracks
    track_indices= position.linear(r).cropped_indices;
    subplot(1,3,1)
    plot(position.clean.x(track_indices),position.clean.y(track_indices),parameters.plot_color_dot{r})
    subplot(1,3,2)
    plot(position.clean.t(track_indices),position.clean.x(track_indices),parameters.plot_color_dot{r})
    subplot(1,3,3)
    plot(position.clean.t(track_indices),position.clean.y(track_indices),parameters.plot_color_dot{r})
end
subplot(1,3,1)
plot(position.clean.x(position.sleep.index),position.clean.y(position.sleep.index),'y.')
axis equal
subplot(1,3,2)
plot(position.clean.t(position.sleep.index),position.clean.x(position.sleep.index),'y.')
subplot(1,3,3)
plot(position.clean.t(position.sleep.index),position.clean.y(position.sleep.index),'y.')

save('position_data','position');
end

%% Plots X and Y position as histogram

function histoplot(y,x) %x and y flipped on purpose
x_bin_width = 1;
y_bin_width = 1;
x_bins = (min(x)-x_bin_width/2):x_bin_width:(max(x)+x_bin_width/2);
y_bins = (min(y)-y_bin_width/2):y_bin_width:(max(y)+x_bin_width/2);
x_bin_centres = min(x):x_bin_width:max(x);
y_bin_centres = min(y):y_bin_width:max(y);
n=histcounts2(x,y,x_bins,y_bins); %calculates  the  amount  of  times  there's  a  repeated X  or  Y  position.
n=log2(n+1); %log  scales  of  histcount - where higher amount of (x or y) repetitions will  be  logged  to (given)  higher  values  (e.g.  log(2 points)=0.3 - but log(200 points)=2.3).
imagesc(y_bin_centres,x_bin_centres,n)
reverse_bone=flipud(colormap(bone));
reverse_bone(:,3)=1;
colormap(reverse_bone);
axis xy
axis equal
end

