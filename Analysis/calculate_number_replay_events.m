
function num_events = calculate_number_replay_events

cd('X:\BendorLab\Drobo\Neural and Behavioural Data\Rate remapping\Data')
load('folders_to_process_remapping.mat')

for f = 1 : length(folders)
    
    cd(folders{f})
    load('sorted_replay_wcorr.mat')
    
    num_events.T1(f) = length(sorted_replay(1).index.track(1).behaviour)+length(sorted_replay(1).index.track(2).behaviour);
    num_events.T2(f) = length(sorted_replay(2).index.track(1).behaviour)+length(sorted_replay(2).index.track(2).behaviour);
    num_events.pre(f) = length(sorted_replay(1).index.PRE)+length(sorted_replay(2).index.PRE);
    num_events.post(f) = length(sorted_replay(1).index.POST)+length(sorted_replay(2).index.POST);
    num_events.T1_local(f) = length(sorted_replay(1).index.track(1).behaviour);
    num_events.T2_local(f) = length(sorted_replay(2).index.track(2).behaviour);
    num_events.pre_T1(f) = length(sorted_replay(1).index.PRE);
    num_events.post_T1(f) = length(sorted_replay(1).index.POST);
    num_events.pre_T2(f) = length(sorted_replay(2).index.PRE);
    num_events.post_T2(f) = length(sorted_replay(2).index.POST);
    cd ..
end

all_events = sum([num_events.T1 num_events.T2 num_events.pre num_events.post]);
num_events.total_events = all_events;

total_events_session = num_events.T1+num_events.T2+num_events.pre+num_events.post;
num_events.mean_total_session = mean(total_events_session);
num_events.std_total_sesion = std(total_events_session);

num_events.T1_mean = mean(num_events.T1);
num_events.T1_std = std(num_events.T1);
num_events.T2_mean = mean(num_events.T2);
num_events.T2_std = std(num_events.T2);


num_events.T1_local_mean = mean(num_events.T1_local);
num_events.T1_local_std = std(num_events.T1_local);
num_events.T2_local_mean = mean(num_events.T2_local);
num_events.T2_local_std = std(num_events.T2_local);

num_events.pre_mean = mean(num_events.pre);
num_events.pre_std = std(num_events.pre);
num_events.post_mean = mean(num_events.post);
num_events.post_std = std(num_events.post);

num_events.pre_T1_mean = mean(num_events.pre_T1);
num_events.pre_T1_std = std(num_events.pre_T1);
num_events.post_T1_mean = mean(num_events.post_T1);
num_events.post_T1_std = std(num_events.post_T1);

num_events.pre_T2_mean = mean(num_events.pre_T2);
num_events.pre_T2_std = std(num_events.pre_T2);
num_events.post_T2_mean = mean(num_events.post_T2);
num_events.post_T2_std = std(num_events.post_T2);


end