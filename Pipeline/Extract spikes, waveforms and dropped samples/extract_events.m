function extract_events
%extracts timestamps of behavioural events marked with a ttl

[data,timestamps,~]=read_neuralynx_file('Events.nev');
events_data.timestamps=timestamps';
events_data.ttl=data.ttl';
events_data.EventStrings=data.EventStrings;
events_data.event_id=data.Event_ID';
events_data.extras=data.Extras;

if ~exist('events_data','var')
    disp('Could not extract events data. Missing Events.nev file inside Raw data')
end

save('extracted_events.mat','events_data','-v7.3');
end