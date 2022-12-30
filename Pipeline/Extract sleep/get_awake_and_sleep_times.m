function [states,state_name]=get_awake_and_sleep_times
load('extracted_position.mat');
load('extracted_sleep_state.mat');
if sleep_state.sleep_indices.start(1) ~= 1 %if the first index is not 1, it means that session starts with rat awake
    awake_indices= [1, sleep_state.sleep_indices.start(1)-1];
    awake_indices= [awake_indices ; (sleep_state.sleep_indices.stop(1:end-1)+1)' , (sleep_state.sleep_indices.start(2:end)-1)'];
else
    awake_indices= [sleep_state.sleep_indices.stop(1:end-1)+1' , sleep_state.sleep_indices.start(2:end)-1'];
end
if sleep_state.sleep_indices.stop(end) ~= length(position.t) % if the rat is not sleeping at the end of the session
    awake_indices= [awake_indices; sleep_state.sleep_indices.stop(end)+1 , length(position.t)];
end
sleep_indices= [sleep_state.sleep_indices.start' , sleep_state.sleep_indices.stop'];
states= [{awake_indices}, {sleep_indices}];
state_name=[{'awake'}, {'sleep'}];
end