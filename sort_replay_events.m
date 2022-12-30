function [sorted_replay,time_range]=sort_replay_events(option,varargin)
% option is [] if no reexposure, and a vector for reexposures, e.g. [ 1 3 ; 2 4 ]
% varargin can be 'wcorr', 'spearman' ,'control_fixed_spike' or 'rate_detection_control'

if isempty(varargin)
    load significant_replay_events;
else
    switch varargin{1}
        case 'wcorr'
            if exist('significant_replay_events_wcorr.mat')==2
                load significant_replay_events_wcorr
%                 significant_replay_events= significant_events_wcorr;
            elseif exist('significant_replay_events_wcorr_individual_exposures.mat')==2
                load significant_replay_events_wcorr_individual_exposures
            else
                disp('wcorr not found, loading default')
                load significant_replay_events;
            end
        case 'spearman'
            if exist('significant_replay_events_spearman.mat')==2
                load significant_replay_events_spearman
%                 significant_replay_events= significant_events_spearman;
            elseif exist('significant_replay_events_spearman_individual_exposures.mat')==2
                load significant_replay_events_spearman_individual_exposures
            else
                disp('spearman not found, loading default')
                load significant_replay_events;
            end
        case 'control_fixed_spike'
            if exist('significant_replay_events_wcorr_FIXED.mat')==2
                load significant_replay_events_wcorr_FIXED
            elseif exist('significant_replay_events_wcorr_individual_exposures_FIXED.mat')==2
                load significant_replay_events_wcorr_individual_exposures_FIXED
            else
                disp('control_fixed_spikes not found, loading default')
                load significant_replay_events_wcorr;
            end
        case 'rate_detection_control'
            if exist('significant_replay_events_wcorr_RATE.mat')==2
                load significant_replay_events_wcorr_RATE
            elseif exist('significant_replay_events_wcorr_individual_exposures_RATE.mat')==2
                load significant_replay_events_wcorr_individual_exposures_RATE
            else
                disp('control_fixed_spikes not found, loading default')
                load significant_replay_events_wcorr;
            end
        otherwise
            disp('did not recognise input, loading default');
            load significant_replay_events;
    end
end

load extracted_position
load extracted_sleep_state


number_of_tracks=length(position.linear);
for track=1:number_of_tracks
    time_range.track(track).behaviour=[min(position.linear(track).timestamps) max(position.linear(track).timestamps)]';
end

time_range.pre=[min(position.t) time_range.track(1).behaviour(1)];
if isempty(option)  %NO REXPOSURE
    time_range.post=[time_range.track(2).behaviour(2) max(position.t)];
else %assumes two tracks
    time_range.post=[time_range.track(2).behaviour(2) time_range.track(3).behaviour(1)];  %when track 2 is finished till the rexposure of track 1 (labelled as 3)
end

sleep_start=position.t(sleep_state.sleep_indices.start);
sleep_stop=position.t(sleep_state.sleep_indices.stop);

time_range.sleep=[sleep_start;sleep_stop];
time_range.awake=[time_range.pre(1) sleep_stop; sleep_start max(position.t)];
total_time_asleep=sum(time_range.sleep(2,:)-time_range.sleep(1,:))/60
total_time_awake=sum(time_range.awake(2,:)-time_range.awake(1,:))/60

% Find sleep and awake times for PRE
index=find(time_range.sleep(1,:)>=time_range.pre(1) & time_range.sleep(1,:)<time_range.pre(2));
time_range.sleepPRE=time_range.sleep(:,index);
index=find(time_range.awake(1,:)>=time_range.pre(1) & time_range.awake(1,:)<time_range.pre(2));
time_range.awakePRE=time_range.awake(:,index);

%adjust for final time when animal is put on track
time_range.awakePRE(find(time_range.awakePRE>time_range.pre(2)))=time_range.pre(2);
time_range.sleepPRE(find(time_range.sleepPRE>time_range.pre(2)))=time_range.pre(2);

% Find sleep and awake times for POST 
index=find(time_range.sleep(2,:)>=time_range.post(1) & time_range.sleep(1,:)<=time_range.post(2));
time_range.sleepPOST=time_range.sleep(:,index);
index=find(time_range.awake(2,:)>=time_range.post(1) & time_range.awake(1,:)<=time_range.post(2));
time_range.awakePOST=time_range.awake(:,index);

%add first awake time when animal is taken off track
time_range.awakePOST(find(time_range.awakePOST<time_range.post(1)))=time_range.post(1);
time_range.sleepPOST(find(time_range.sleepPOST<time_range.post(1)))=time_range.post(1);
time_range.awakePOST(find(time_range.awakePOST>time_range.post(2)))=time_range.post(2);
time_range.sleepPOST(find(time_range.sleepPOST>time_range.post(2)))=time_range.post(2);

time_range.PRE_CUMULATIVE=compute_cumulative_time(time_range.pre');
time_range.awakePRE_CUMULATIVE=compute_cumulative_time(time_range.awakePRE);
time_range.sleepPRE_CUMULATIVE=compute_cumulative_time(time_range.sleepPRE);
time_range.POST_CUMULATIVE=compute_cumulative_time(time_range.post');
time_range.awakePOST_CUMULATIVE=compute_cumulative_time(time_range.awakePOST);
time_range.sleepPOST_CUMULATIVE=compute_cumulative_time(time_range.sleepPOST);
for track=1:number_of_tracks
    time_range.RUN_CUMULATIVE{track} = compute_cumulative_time(time_range.track(track).behaviour);
end

for track=1:number_of_tracks
    for j=1:number_of_tracks
        sorted_replay(track).index.track(j).behaviour=[];
        sorted_replay(track).ref_index.track(j).behaviour=[];
    end
    sorted_replay(track).index.PRE=[];
    sorted_replay(track).ref_index.PRE=[];
    sorted_replay(track).index.awakePRE=[];
    sorted_replay(track).ref_index.awakePRE=[];
    sorted_replay(track).index.sleepPRE=[];
    sorted_replay(track).ref_index.sleepPRE=[];
    sorted_replay(track).index.POST=[];
    sorted_replay(track).ref_index.POST=[];
    sorted_replay(track).index.awakePOST=[];
    sorted_replay(track).ref_index.awakePOST=[];
    sorted_replay(track).index.sleepPOST=[];
    sorted_replay(track).ref_index.sleepPOST=[];
    for j=1:number_of_tracks
        sorted_replay(track).event_time.track(j).behaviour=[];
    end
    sorted_replay(track).event_time.PRE=[];
    sorted_replay(track).event_time.awakePRE=[];
    sorted_replay(track).event_time.sleepPRE=[];
    sorted_replay(track).event_time.POST=[];
    sorted_replay(track).event_time.awakePOST=[];
    sorted_replay(track).event_time.sleepPOST=[];
end

for track=1:number_of_tracks
    event_times=significant_replay_events.track(track).event_times;
    for j=1:length(event_times)
        index = j;
        ref_index = significant_replay_events.track(track).ref_index(j);
        for k=1:number_of_tracks
            if check_if_event_is_in_time_window(event_times(j),time_range.track(k).behaviour)
                sorted_replay(track).index.track(k).behaviour=[sorted_replay(track).index.track(k).behaviour index];
                sorted_replay(track).ref_index.track(k).behaviour=[sorted_replay(track).ref_index.track(k).behaviour ref_index];
                sorted_replay(track).event_time.track(k).behaviour=[sorted_replay(track).event_time.track(k).behaviour event_times(j)];
            end
        end
        if check_if_event_is_in_time_window(event_times(j),time_range.pre')
            sorted_replay(track).index.PRE=[sorted_replay(track).index.PRE index];
            sorted_replay(track).ref_index.PRE=[sorted_replay(track).ref_index.PRE ref_index];
            sorted_replay(track).event_time.PRE=[sorted_replay(track).event_time.PRE event_times(j)];
        end
        if check_if_event_is_in_time_window(event_times(j),time_range.awakePRE)
            sorted_replay(track).index.awakePRE=[sorted_replay(track).index.awakePRE index];
            sorted_replay(track).ref_index.awakePRE=[sorted_replay(track).ref_index.awakePRE ref_index];
            sorted_replay(track).event_time.awakePRE=[sorted_replay(track).event_time.awakePRE event_times(j)];
        end
        if check_if_event_is_in_time_window(event_times(j),time_range.sleepPRE)
            sorted_replay(track).index.sleepPRE=[sorted_replay(track).index.sleepPRE index];
            sorted_replay(track).ref_index.sleepPRE=[sorted_replay(track).ref_index.sleepPRE ref_index];
            sorted_replay(track).event_time.sleepPRE=[sorted_replay(track).event_time.sleepPRE event_times(j)];
        end
        if check_if_event_is_in_time_window(event_times(j),time_range.post')
            sorted_replay(track).index.POST=[sorted_replay(track).index.POST index];
            sorted_replay(track).ref_index.POST=[sorted_replay(track).ref_index.POST ref_index];
            sorted_replay(track).event_time.POST=[sorted_replay(track).event_time.POST event_times(j)];
        end
        if check_if_event_is_in_time_window(event_times(j),time_range.awakePOST)
            sorted_replay(track).index.awakePOST=[sorted_replay(track).index.awakePOST index];
            sorted_replay(track).ref_index.awakePOST=[sorted_replay(track).ref_index.awakePOST ref_index];
            sorted_replay(track).event_time.awakePOST=[sorted_replay(track).event_time.awakePOST event_times(j)];
        end
        if check_if_event_is_in_time_window(event_times(j),time_range.sleepPOST)
            sorted_replay(track).ref_index.sleepPOST=[sorted_replay(track).ref_index.sleepPOST ref_index];
            sorted_replay(track).index.sleepPOST=[sorted_replay(track).index.sleepPOST index];
            sorted_replay(track).event_time.sleepPOST=[sorted_replay(track).event_time.sleepPOST event_times(j)];
        end
    end
end

for track=1:number_of_tracks
    sorted_replay(track).cumulative_event_time.PRE=interpolate_cumulative_time(time_range.PRE_CUMULATIVE,...
        time_range.pre',sorted_replay(track).event_time.PRE);
    sorted_replay(track).cumulative_event_time.awakePRE=interpolate_cumulative_time(time_range.awakePRE_CUMULATIVE,...
        time_range.awakePRE,sorted_replay(track).event_time.awakePRE);
    sorted_replay(track).cumulative_event_time.sleepPRE=interpolate_cumulative_time(time_range.sleepPRE_CUMULATIVE,...
        time_range.sleepPRE,sorted_replay(track).event_time.sleepPRE);
    sorted_replay(track).cumulative_event_time.POST=interpolate_cumulative_time(time_range.POST_CUMULATIVE,...
        time_range.post',sorted_replay(track).event_time.POST);
    sorted_replay(track).cumulative_event_time.awakePOST=interpolate_cumulative_time(time_range.awakePOST_CUMULATIVE,...
        time_range.awakePOST,sorted_replay(track).event_time.awakePOST);
    sorted_replay(track).cumulative_event_time.sleepPOST=interpolate_cumulative_time(time_range.sleepPOST_CUMULATIVE,...
        time_range.sleepPOST,sorted_replay(track).event_time.sleepPOST);
    for k=1:number_of_tracks
        sorted_replay(track).cumulative_event_time.track(k).behaviour=interpolate_cumulative_time(time_range.RUN_CUMULATIVE{k},...
            time_range.track(k).behaviour,sorted_replay(track).event_time.track(k).behaviour);
    end
end

[sorted_replay(:).method]= deal(varargin{1});

save time_range time_range
if isempty(varargin)
    save sorted_replay sorted_replay
else
    switch varargin{1}
        case 'wcorr'
            save sorted_replay_wcorr sorted_replay
        case 'spearman'
            save sorted_replay_spearman sorted_replay
        case 'control_fixed_spike'
            save sorted_replay_wcorr_FIXED sorted_replay
        case 'rate_detection_control'
            save sorted_replay_wcorr_RATE sorted_replay
        otherwise
            save sorted_replay sorted_replay
    end
end

end

function output=check_if_event_is_in_time_window(event_time,time_window)
output=0;
for i=1:size(time_window,2)
    if event_time>=time_window(1,i) & event_time<time_window(2,i)
        output=1;
    end
end
end

function cumulative_time=compute_cumulative_time(time)
% Takes start and end times for a specific period (e.g. PRE sleep), normalized it to the start of the sleep epoch within the sleep period and 
% then calculates cummulative sum for the start and end timestamps (stacks all epochs together)

cumulative_time=time-time(1,:); %normalize by epoch start
cumulative_time(2,:)=cumsum(cumulative_time(2,:));%cummulative sum of stop times
cumulative_time(1,2:end)=cumulative_time(2,1:(end-1));
end

function cumulative_event_times=interpolate_cumulative_time(cumulative_time,time, event_times)
% Interpolate replay event times to cumulative time

if isempty(time)
    cumulative_event_times=NaN;
else
    time(1,2:end)=time(1,2:end)+1e-10;
    cumulative_time(1,2:end)=cumulative_time(1,2:end)+1e-10;
    t=NaN(1,size(time,1)*size(time,2));
    t_c=NaN(1,size(time,1)*size(time,2));
    t(1:2:end)=time(1,:);
    t(2:2:end)=time(2,:);
    t_c(1:2:end)=cumulative_time(1,:);
    t_c(2:2:end)=cumulative_time(2,:);
    cumulative_event_times=interp1(t,t_c,event_times,'linear');
end
end
