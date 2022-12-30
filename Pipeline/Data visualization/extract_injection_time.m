function injection_time=extract_injection_time


[data,timestamps,info]=read_neuralynx_file('Events.nev');
if ~isempty(data)
%FOR DREADDS%%%%%%%%%%%%%%%%%%%%%%%%%%%%
injection_time=1e-6*CNO_injection(timestamps,data.EventStrings);
else
    injection_time=[];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(injection_time)
    load extracted_position
   injection_time=min(position.t); 
end
end



function events=extract_events
[Timestamps, EventIDs, TTLs, Extras, EventStrings, Header] = Nlx2MatEV('Events.nev', [1 1 1 1 1],1,1,[]);
events.timestamps=Timestamps';
events.ttl=TTLs';
events.EventStrings=EventStrings;
end

function injection_time=CNO_injection(timestamps,EventStrings)
injection_time=[];
for i=1:length(EventStrings)
    if (~isempty(strfind(EventStrings{i},'CNO')) | ~isempty(strfind(EventStrings{i},'PBS')) | ...
            ~isempty(strfind(EventStrings{i},'cno')) | ~isempty(strfind(EventStrings{i},'pbs')) | ...
            ~isempty(strfind(EventStrings{i},'VEHICLE')) | ~isempty(strfind(EventStrings{i},' saline')));
        injection_time=timestamps(i);
        break;
    end
end
end