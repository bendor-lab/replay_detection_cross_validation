% EXTRACTS SLEEP STATE
% Considers that animal is sleeping when:
% Velocity is less than 4cm/s  within a 60 second window
% Median firing rate across neurons with a lower firing rate (lower 1/3 of population) has a zscore > 0 within a 60 second window
% Interpolates replay event time with time_bin_centre (using 'nearest' option) to determine if the animal is sleep (sleep = 1) or not (sleep = -1)
% Loads extracted_clusters.m & extracted_position.m
% varargin input allows mua and speed threshold to be manually adjusted

function sleep_state = sleep_stager(function_mode)

load extracted_clusters
load extracted_position
parameters = list_of_parameters;
response='y';   %stay in while loop till this is changed to 'n'
time_bin_edges = min(clusters.spike_times):60:max(clusters.spike_times);  % 60 second bins
time_bin_centres = time_bin_edges(1:end-1)+30;
v_cm = smooth(position.v_cm,25*60+1);  %60 seconds smoothing
velocity = interp1(position.t,v_cm,time_bin_centres,'linear'); % interpolate velocity to bins timestamps

% For each unit, create spike count for each time bin and normalize it
unique_cluster_spikes_id=unique(clusters.spike_id);
for i = 1 : length(unique_cluster_spikes_id)
    spikes(i,:) = histcounts(clusters.spike_times(clusters.spike_id==unique_cluster_spikes_id(i)),time_bin_edges);
    spikes_sum(i) = sum(spikes(i,:)); %sum of all the spikes for that unit
    spikes(i,:) = zscore(spikes(i,:));
end

% Sort units by number of spikes
[~,index] = sort(spikes_sum);
units_with_more_spikes = spikes(index(1:ceil(length(index)/3)),:); % find top 1/3 of units with more spikes
units_with_more_spikes(units_with_more_spikes>2) = 2;  %cap max and min values
units_with_more_spikes(units_with_more_spikes<-2) = -2;


while(strcmp(response,'y'))
    %default values:
    MUA_THRESHOLD = 0;
    SPEED_THRESHOLD = parameters.speed_threshold;
    % Raster plot from top third of units with more spikes
    figure
    ax(1)=subplot(4,1,1:2);
    imagesc(time_bin_edges,1:size(units_with_more_spikes,1),units_with_more_spikes)
    axis xy
    colormap(jet)
    ylabel('most active units')
    ax(2)=subplot(4,1,3);
    plot(time_bin_centres,zscore(median(units_with_more_spikes,1)),'b'); %zscored of number of spikes per time bin (median of all units)
    hold on
    velocity(velocity>5) = 5; %cap velocity for plotting
    plot(time_bin_centres,(velocity),'r');
    plot([min(time_bin_centres) max(time_bin_centres)],[MUA_THRESHOLD MUA_THRESHOLD],'Color',[0, 0.4470, 0.7410]) %MUA threshold
    plot([min(time_bin_centres) max(time_bin_centres)],[SPEED_THRESHOLD SPEED_THRESHOLD],'Color',[0.6350, 0.0780, 0.1840]) %velocity threshold
    
    ylabel('zscored spikes / zscored speed')
    legend('zscored number of spikes per time bin','velocity')
    
    if strcmp(function_mode,'manual') %if user wants to manually change thresholds
    new_MUA_thresh = input('If needed, input new MUA_threshold(default [0])','s');
    new_speed_thresh = input('If needed, input new Speed_threshold (default [4])','s');
    if ~isempty(new_MUA_thresh)
        MUA_THRESHOLD = str2double(new_MUA_thresh);
    end
    if ~isempty(new_speed_thresh)
        SPEED_THRESHOLD = str2double(new_speed_thresh);
    end
    end
    
    plot([min(time_bin_centres) max(time_bin_centres)],[MUA_THRESHOLD MUA_THRESHOLD],'Color',[0, 0.4470, 0.7410],'LineStyle','- -') %MUA threshold
    plot([min(time_bin_centres) max(time_bin_centres)],[SPEED_THRESHOLD SPEED_THRESHOLD],'Color',[0.6350, 0.0780, 0.1840],'LineStyle','- -') %velocity threshold
    
    ax(3)=subplot(4,1,4);
    sleep = -ones(size(median(units_with_more_spikes,1)));
    sleep(find(zscore(median(units_with_more_spikes,1)) > MUA_THRESHOLD & velocity < SPEED_THRESHOLD)) = 1; % find when there's spikes and low velocity
    bar(time_bin_centres,sleep);
    xlabel('time')
    ylabel(' awake   -   sleep')
    linkaxes(ax,'x');
    
    % Save in structure
    sleep_state.time = time_bin_centres;
    sleep_state.state_binned = sleep; % 1 for sleep, -1 for awake
    
    % Find start and end of sleep periods as position.t indices
    states = interp1(sleep_state.time,sleep_state.state_binned,position.t,'nearest'); %interpolate back to position.t
    states(isnan(states))= -1; % Replace NaNs at start and end for -1 (non sleep)
    sleep_state.state = states;
    
    state_change_indxs = find(abs(diff(states)) > 1);
    sleep_state.sleep_indices.start= state_change_indxs(1:2:end-1);
    sleep_state.sleep_indices.stop= state_change_indxs(2:2:end);
    
    % Calculate time spent sleeping based on velocity for PRE/POTS/POST
    pre_idx = sleep_state.time < min(position.linear(1).timestamps);
    pre_sleep_idx = find(sleep_state.state_binned(pre_idx) > 0);
    sleep_state.time_slept.PRE = length(pre_sleep_idx)*60; % in seconds
    
    for i = 1:length(position.linear)-1
        pot_idx = sleep_state.time > max(position.linear(i).timestamps) & sleep_state.time < min(position.linear(i+1).timestamps);
        pot_sleep_idx = find(sleep_state.state_binned(pot_idx) > 0);
        if length(position.linear)-1 == 1 || length(position.linear)-1 == 2 % If 2 or 3 track exposures in the session
            sleep_state.time_slept.(['sleep_pot' num2str(i)]) = length(pot_sleep_idx)*60; % in seconds
        elseif length(position.linear)-1 == 3 % If 4 track exposures in the session
            if i == 1 || i == 3
                sleep_state.time_slept.(['sleep_pot' num2str(i)]) = length(pot_sleep_idx)*60; % in seconds
            elseif i == 2
                sleep_state.time_slept.(['INTER_POST' num2str(i)]) = length(pot_sleep_idx)*60; % in seconds
            end
        end
    end
    
    post_idx = sleep_state.time > max(position.linear(length(position.linear)).timestamps);
    post_sleep_idx = find(sleep_state.state_binned(post_idx) > 0);
    sleep_state.time_slept.FINAL_POST = length(post_sleep_idx)*60; % in seconds
    sleep_state.legend = {'Sleep_state', 1, 'Awake_state', -1};
    sleep_state.MUA_THRESHOLD = MUA_THRESHOLD;
    sleep_state.SPEED_THRESHOLD = SPEED_THRESHOLD;
    
    
    if strcmp(function_mode,'manual') %if user wants to manually change thresholds
    response = input('Do you want to refine your thresholds: y or n:   ','s');
    else
        response='n';
    end
end
save extracted_sleep_state sleep_state

end
