
function replay_decoding_split_events
%splits each replay event approximately in half, based on when the MUA is
%lowest in the middle third of the event.  This is used to seperately analyze the
%first and second half of long duration replay events

% Load parameters
parameters = list_of_parameters;
load('decoded_replay_events')
load('extracted_replay_events')

num_replay_events = length(replay.onset);

for track = 1:length(decoded_replay_events)
    for event = 1 : num_replay_events
        spikes=decoded_replay_events(track).replay_events(event).spikes(:,2);  %spike times
         %find time bin where midpoint of replay event occurs
         [~,index]=min(abs(decoded_replay_events(track).replay_events(event).timebins_centre-replay.midpoint(event))); 
        bin_after_index=index+1; % to allow segments to overlap by one bin
        if bin_after_index>length(decoded_replay_events(track).replay_events(event).timebins_edges)
            bin_after_index=length(decoded_replay_events(track).replay_events(event).timebins_edges);
        end
        
        %each event split into two segments: decoded_replay_events1 and decoded_replay_events2
        %first segment
        decoded_replay_events1(track).replay_events(event).replay_id=decoded_replay_events(track).replay_events(event).replay_id;
        decoded_replay_events1(track).replay_events(event).midpoint=replay.midpoint(event);
        decoded_replay_events1(track).replay_events(event).onset=replay.onset(event);
        decoded_replay_events1(track).replay_events(event).offset=replay.midpoint(event);
        decoded_replay_events1(track).replay_events(event).duration=replay.midpoint(event)-replay.onset(event);
        spike_index1=find(spikes<=replay.midpoint(event));
        decoded_replay_events1(track).replay_events(event).spikes(:,1)=decoded_replay_events(track).replay_events(event).spikes(spike_index1,1);
        decoded_replay_events1(track).replay_events(event).spikes(:,2)=decoded_replay_events(track).replay_events(event).spikes(spike_index1,2);
         decoded_replay_events1(track).replay_events(event).timebins_edges=decoded_replay_events(track).replay_events(event).timebins_edges(1:bin_after_index);
        decoded_replay_events1(track).replay_events(event).timebins_centre=decoded_replay_events(track).replay_events(event).timebins_centre(1:index);
        decoded_replay_events1(track).replay_events(event).timebins_index=1:index;
        decoded_replay_events1(track).replay_events(event).decoded_position=decoded_replay_events(track).replay_events(event).decoded_position(:,1:index);
        
        %second segment
        decoded_replay_events2(track).replay_events(event).replay_id=decoded_replay_events(track).replay_events(event).replay_id;
        decoded_replay_events2(track).replay_events(event).midpoint=replay.midpoint(event);
        decoded_replay_events2(track).replay_events(event).onset=replay.midpoint(event);
        decoded_replay_events2(track).replay_events(event).offset=replay.offset(event);
        decoded_replay_events2(track).replay_events(event).duration=replay.offset(event)-replay.midpoint(event);
        spike_index2=find(spikes>=replay.midpoint(event));
        decoded_replay_events2(track).replay_events(event).spikes(:,1)=decoded_replay_events(track).replay_events(event).spikes(spike_index2,1);
        decoded_replay_events2(track).replay_events(event).spikes(:,2)=decoded_replay_events(track).replay_events(event).spikes(spike_index2,2);
        decoded_replay_events2(track).replay_events(event).timebins_edges=decoded_replay_events(track).replay_events(event).timebins_edges(index:end);
        decoded_replay_events2(track).replay_events(event).timebins_centre=decoded_replay_events(track).replay_events(event).timebins_centre(index:end);
        decoded_replay_events2(track).replay_events(event).timebins_index=index:length(decoded_replay_events(track).replay_events(event).timebins_centre);
        decoded_replay_events2(track).replay_events(event).decoded_position=decoded_replay_events(track).replay_events(event).decoded_position(:,index:end);
    end
end
% Saves structure
save decoded_replay_events_segments decoded_replay_events1 decoded_replay_events2
end


