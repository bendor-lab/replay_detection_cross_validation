function [] = proportion_of_unique_place_fields(folders,option,method)

load('Tables\subsets_of_cells.mat')

total_number{1}(1:3) = 0;
total_number{2}(1:3) = 0;

for nfolder = 1:length(folders) % Get total number of replay events
    for nshuffle = 1:2 % 1 is original and 2 is global remapped
        session_total_number{nshuffle}(nfolder,1:3) = 0;
    end
end
unique_good_cell_index = [];
common_good_cell_index = [];

for nfolder = 1:length(folders) % Get total number of replay events
    cd(folders{nfolder})
    load extracted_place_fields_BAYESIAN.mat
    for nshuffle = 1:2 % 1 is original and 2 is global remapped


        common_good_cell_index{nfolder} = subset_of_cells.cell_IDs{1}...
            (intersect(find(subset_of_cells.cell_IDs{1} >= 1000*nfolder),find(subset_of_cells.cell_IDs{1} <= 1000*(nfolder+1))))- nfolder*1000;

        unique_good_cell_index{nfolder} = place_fields_BAYESIAN.good_place_cells(~ismember(place_fields_BAYESIAN.good_place_cells,common_good_cell_index{nfolder}));

        if nshuffle == 2
            cd global_remapped_shuffles
            number_of_global_remapped_shuffles = length(dir('shuffle_*'));
            cd ..
        elseif nshuffle == 1 % Original data, no shuffle
            number_of_global_remapped_shuffles = 1;
        end

        load decoded_replay_events_segments

        [~,time_range]=sort_replay_events([],'wcorr');

        for event = 1:length(decoded_replay_events1(1).replay_events)

            event_time = decoded_replay_events1(1).replay_events(event).midpoint;

            if event_time >= time_range.pre(1) & event_time <= time_range.pre(2) %If PRE

                total_number{nshuffle}(1) = total_number{nshuffle}(1) + 1*number_of_global_remapped_shuffles;
                session_total_number{nshuffle}(nfolder,1) = session_total_number{nshuffle}(nfolder,1) + 1*number_of_global_remapped_shuffles;

            elseif event_time >= time_range.post(1) & event_time <= time_range.post(2) %If POST

                total_number{nshuffle}(3) = total_number{nshuffle}(3) + 1*number_of_global_remapped_shuffles;
                session_total_number{nshuffle}(nfolder,3) = session_total_number{nshuffle}(nfolder,3) + 1*number_of_global_remapped_shuffles;

            elseif event_time >= time_range.track(1).behaviour(1) & event_time <= time_range.track(1).behaviour(2)
                total_number{nshuffle}(2) = total_number{nshuffle}(2) + 1*number_of_global_remapped_shuffles;
                session_total_number{nshuffle}(nfolder,2) = session_total_number{nshuffle}(nfolder,2) + 1*number_of_global_remapped_shuffles;

            elseif event_time >= time_range.track(2).behaviour(1) & event_time <= time_range.track(2).behaviour(2)
                total_number{nshuffle}(2) = total_number{nshuffle}(2) + 1*number_of_global_remapped_shuffles;
                session_total_number{nshuffle}(nfolder,2) = session_total_number{nshuffle}(nfolder,2) + 1*number_of_global_remapped_shuffles;

            end
        end

    end
    cd ..
end

load('P:\ground_truth_replay_analysis\Dropbo_data8\ground_truth_original\log_odd_difference_multiple_shuffles.mat')
index = [];
epoch_index = [];
cd ground_truth_original
load log_odd_wcorr_PRE_place_POST_time
cd ..


% Get index for PRE,  RUN, POST
% states = [-1 0 1 2 3 4 5]; % sleepPRE, awakePRE, RUN T1, RUN T2, sleepPOST and awakePOST and multi-track
% states = [-1 0 1 2 3 4]; % sleepPRE, awakePRE, RUN T1, RUN T2, sleepPOST and awakePOST
states = [-1 0 1 2];

for k=1:length(states) % sort track id and behavioral states (pooled from 5 sessions)
    % Find intersect of behavioural state and ripple peak threshold
    state_index = find(log_odd.behavioural_state==states(k));

    if k == 1 % PRE
        epoch_index{1} = state_index;
    elseif k == 2 % POST
        epoch_index{3} = state_index;
    elseif k == 3 % RUN Track 1
        epoch_index{2} = state_index;
    elseif k == 4 % RUN Track 2
        epoch_index{2} = [epoch_index{2} state_index];
    end
end


%     colour_line2 = {'k--','b--','r--','g--'};
colour_symbol={'bo','ro','go','ko'};
Behavioural_epoches = {'PRE','RUN','POST'};
end

cell_id_event = [];
count = 1;
for nfolder = 1:length(folders)
    cd(folders{nfolder})
    load decoded_replay_events
    load decoded_replay_events_segments
    load extracted_replay_events
    for event = 1:length(decoded_replay_events(1).replay_events)
        if replay.ripple_peak(event) > 3
            cell_id_event{count}{1} = unique(decoded_replay_events(1).replay_events(event).spikes(:,1));
            cell_id_event{count}{2} = unique(decoded_replay_events1(1).replay_events(event).spikes(:,1));
            cell_id_event{count}{3} = unique(decoded_replay_events2(1).replay_events(event).spikes(:,1));
            count = count + 1;
        end
    end
    cd ..
end

p_val_threshold = round(10.^[-3:0.01:log10(0.2)],4);
p_val_threshold = flip(unique(p_val_threshold));
p_val_threshold(1) = 0.2;
low_threshold = find(ismembertol(p_val_threshold,[0.0501 0.02 0.01 0.005 0.002 0.001],eps) == 1);

for epoch = 1:3
    for condition = 1:2
        % at p <= 0.05
        if condition == 1
            threshold = 0.05;
        else
            % at false positive rate ~0.05
            [c index(epoch)] = min(abs(mean(percent_shuffle_events{4}{epoch},2) - 0.05));
            matching_threshold(epoch) = p_val_threshold(index(epoch))
            threshold = p_val_threshold(index(epoch));
        end
        nshuffle = 1;

        events = epoch_index{epoch};
        
        
        %     for track = 1:2
        this_segment_event = events(find(log_odd.segment_id(1,events) == 1));
        track_1_index = this_segment_event(find(log_odd.pvalue(1,this_segment_event) <= threshold));
        this_segment_event = events(find(segment_id{nmethod}{nshuffle}(1,events) > 1));
        track_1_index = [track_1_index this_segment_event(find(log_odd.pvalue(1,this_segment_event) <= threshold/2))];
        
        for event = 1:length(track_1_index)
             experiment = log_odd.experiment(track_1_index(event));
             cell_this_event = cell_id_event{track_1_index(event)}{log_odd.segment_id(1,track_1_index(event))};
             common_good_proportion{epoch}{1}(event) = sum(ismember(cell_this_event,common_good_cell_index{experiment}))/length(cell_this_event);
             unique_good_proportion{epoch}{1}(event) = sum(ismember(cell_this_event,unique_good_cell_index{experiment}))/length(cell_this_event);
        end
        

        this_segment_event = events(find(log_odd.segment_id(2,events) == 1));
        track_2_index = this_segment_event(find(log_odd.pvalue(2,this_segment_event) <= threshold));
        this_segment_event = events(find(log_odd.segment_id(2,events) > 1));
        track_2_index = [track_2_index this_segment_event(find(log_odd.pvalue(2,this_segment_event) <= threshold/2))];
        
        for event = 1:length(track_2_index)
             experiment = log_odd.experiment(track_2_index(event));
             cell_this_event = cell_id_event{track_2_index(event)}{log_odd.segment_id(2,track_2_index(event))};
             common_good_proportion{epoch}{2}(event) = sum(ismember(cell_this_event,common_good_cell_index{experiment}))/length(cell_this_event);
             unique_good_proportion{epoch}{2}(event) = sum(ismember(cell_this_event,unique_good_cell_index{experiment}))/length(cell_this_event);
        end

%         if length(track_1_index)>length(track_2_index)
%             track_1_index(~ismember(track_1_index,track_2_index)
%         else
%             ~ismember(track_2_index,track_1_index)
%         end
    end
end

figure
for epoch = 1:3
    subplot(2,2,epoch)
    histogram(common_good_proportion{epoch}{1},20)
    hold on 
    histogram(unique_good_proportion{epoch}{1},20)
    legend('Common','Unique')
end
sgtitle('At p value 0.05')
% Behavioural_epoches = {'PRE','RUN','POST'};
figure
for epoch = 1:3
    subplot(2,2,epoch)
    histogram(common_good_proportion{epoch}{2},20)
    hold on 
    histogram(unique_good_proportion{epoch}{2},20)
    legend('Common','Unique')
    title(Behavioural_epoches{epoch})
end
sgtitle('At matching false positive rate')
