function [log_odd] = extract_ground_truth_info_backup(folders,option,method,shuffles,ripple_zscore_threshold)
% Input:
% Folders: Each Folder is data for a session
% method = 'linear' or 'wcorr' or 'spearman' or 

% Output:
% log_odd: save the log odd - log(sum(probability of T1)/sum(probability of
% T2)) of the events that are considered significant according to the original bayesian decoding
% as well as the 1000 T1&T2 turning curve shuffles for each event
% It will also contain the session ID, 'ground truth' replay track ID and
% behavioral epoch (-1: PRE, 0: POST, 1: T1 RUN, 2: T2 RUN) and etc

c = 1;
current_directory=pwd;
log_odd=[];
load subsets_of_cells;

cd ground_truth_original
if strcmp(option,'original') || strcmp(option,'common')
    load jump_distance
elseif strcmp(option,'global remapped')
    load jump_distance_global_remapped
end
cd ..

for f = 1:length(folders)
    
    cd(folders{f})
    parameters = list_of_parameters;
    place_fields_BAYESIAN = [];
%     load('significant_replay_events_wcorr.mat');
    load('extracted_place_fields_BAYESIAN.mat');
    load('extracted_position.mat')
    load('extracted_clusters.mat');
    load('extracted_replay_events.mat');
    
    % Stable common good cells (based on subsets_of_cells)
    %     place_cell_index = subset_of_cells.cell_IDs{8}...
    %         (intersect(find(subset_of_cells.cell_IDs{8} >= 1000*f),find(subset_of_cells.cell_IDs{8} <= 1000*(f+1))))- f*1000;
    
    if strcmp(option,'original') || strcmp(option,'common')
        load probability_ratio_original
        load probability_ratio_common_good
    elseif strcmp(option,'global remapped')
        load global_remapped_original_probability_ratio
        load global_remapped_common_good_probability_ratio
    end
    
    
    
    %% Ground Truth Information (Bayesian bias, p value, best score, probability ratio)
    
    
    if strcmp(option,'original') || strcmp(option,'common')

        load scored_replay
        load scored_replay_segments
        
        % gather p value and replay score information
        [p_values, replay_scores] = extract_score_and_pvalue_ground_truth(scored_replay, scored_replay1, scored_replay2, method, []);
                 
        [significant_replay_events sig_event_info] = number_of_significant_replays_ground_truth(0.05,ripple_zscore_threshold,method,shuffles,[]);
    
        % gather time information for sorting the event as PRE, POST or RUN
        [~,time_range]=sort_replay_events([],'wcorr');        
        
        for track=1:2
            for i=1:length(significant_replay_events.track(track).p_value)
                if ~ismember(significant_replay_events.track(track).ref_index(i),significant_replay_events.track(1).ref_index_MULTI) % if not multi-track
                    % save replay event index
                    log_odd.index(c) = significant_replay_events.track(track).ref_index(i);
                    log_odd.multi_track(c) = 0;

                    log_odd.pvalue(c) = sig_event_info.p_value(track,significant_replay_events.track(track).index(i));
                    log_odd.segment_id(c) = sig_event_info.segment_id(track,significant_replay_events.track(track).index(i));

                    % Original 20ms
                    log_odd.normal.probability_ratio(c) = probability_ratio_original{1}(track,log_odd.index(c),log_odd.segment_id(c));
                    log_odd.normal.T1_T2_ratio(c) = probability_ratio_original{1}(1,log_odd.index(c),log_odd.segment_id(c));
                    
                    % True/False Ratio
                    for nshuffle = 1:1000
                        log_odd.normal.probability_ratio_shuffled{c}(nshuffle) = probability_ratio_original{2}{nshuffle}(track,log_odd.index(c),log_odd.segment_id(c));
                    end
                    
                    % Calculate and save zscored log odd
                    data = log(log_odd.normal.probability_ratio(c));
                    shuffled_data = log(log_odd.normal.probability_ratio_shuffled{c});
                    tempt = (data-mean(shuffled_data))/std(shuffled_data);
                    log_odd.normal_zscored.original(1,c) = tempt;
                    
                    % T1/T2 ratio
                    for nshuffle = 1:1000
                        log_odd.normal.T1_T2_ratio_shuffled{c}(nshuffle) = probability_ratio_original{2}{nshuffle}(1,log_odd.index(c),log_odd.segment_id(c));
                    end
                    
                    % Calculate and save zscored log odd
                    data = log(log_odd.normal.T1_T2_ratio(c));
                    shuffled_data = log(log_odd.normal.T1_T2_ratio_shuffled{c});
                    tempt = (data-mean(shuffled_data))/std(shuffled_data);
                    log_odd.normal_zscored.original(2,c) = tempt;
                    
                    
                    % Original 20ms
                    log_odd.common.probability_ratio(c) = probability_ratio_common_good{1}(track,log_odd.index(c),log_odd.segment_id(c));
                    log_odd.common.T1_T2_ratio(c) = probability_ratio_common_good{1}(1,log_odd.index(c),log_odd.segment_id(c));
                    
                    % True/False Ratio
                    for nshuffle = 1:1000
                        log_odd.common.probability_ratio_shuffled{c}(nshuffle) = probability_ratio_common_good{2}{nshuffle}(track,log_odd.index(c),log_odd.segment_id(c));
                    end
                    
                    % Calculate and save zscored log odd
                    data = log(log_odd.common.probability_ratio(c));
                    shuffled_data = log(log_odd.common.probability_ratio_shuffled{c});
                    tempt = (data-mean(shuffled_data))/std(shuffled_data);
                    log_odd.common_zscored.original(1,c) = tempt;
                    
                    % T1/T2 ratio
                    for nshuffle = 1:1000
                        log_odd.common.T1_T2_ratio_shuffled{c}(nshuffle) = probability_ratio_common_good{2}{nshuffle}(1,log_odd.index(c),log_odd.segment_id(c));
                    end
                    
                    % Calculate and save zscored log odd
                    data = log(log_odd.common.T1_T2_ratio(c));
                    shuffled_data = log(log_odd.common.T1_T2_ratio_shuffled{c});
                    tempt = (data-mean(shuffled_data))/std(shuffled_data);
                    log_odd.common_zscored.original(2,c) = tempt;


                    %                     log_odd.track1_pvalue(c)=log10(9e-4+  max(p_values.WHOLE(1,log_odd.index(c),:)));
                    %                     log_odd.track1_pvalue(c)=log10(sig_event_info.p_value(1,log_odd.index(c)));
                    %                     log_odd.pvalue(c)=log10(9e-4+ significant_replay_events.track(track).p_value(i));


                    log_odd.best_score(c)=significant_replay_events.track(track).replay_score(i);
%                     log_odd.track1_best_score(c)=replay_scores.WHOLE(1,log_odd.index(c));

                    log_odd.neuron_count(c)= replay.neuron_count(log_odd.index(c));
                    log_odd.spike_count(c)= replay.spike_count(log_odd.index(c));

                    log_odd.event_duration(c) = significant_replay_events.track(track).event_duration(i);
                    log_odd.ripple_peak(c) = replay.ripple_peak(log_odd.index(c));

                    if isempty(jump_distance{f}{track}{log_odd.index(c)})
                        log_odd.max_jump(c) = nan;
                    else
                        log_odd.max_jump(c) = jump_distance{f}{track}{log_odd.index(c)};
                    end


                    event_time = (replay.onset(log_odd.index(c)) + replay.offset(log_odd.index(c)))/2; % Median event Time
                 
%                     event_time = significant_replay_events.track(track).event_times(i);
                    
                    if event_time <= time_range.pre(2) %If PRE
                        
                        log_odd.behavioural_state(c)= -1;
                        
                    elseif event_time>=time_range.post(1) %If POST
                        
                        log_odd.behavioural_state(c)= 0;
                        
                        
                    elseif event_time>=time_range.track(1).behaviour(1) & event_time<=time_range.track(1).behaviour(2)
                        log_odd.behavioural_state(c)=1;  % If Run T1
                    elseif event_time>=time_range.track(2).behaviour(1) & event_time<=time_range.track(2).behaviour(2)
                        log_odd.behavioural_state(c)=2;  % If Run T2
                        
                        
                        %
                        %                 if event_time < time_range.pre(2) % If Pre
                        %                     for k = 1:size(time_range.sleepPRE,2)
                        %                         if event_time < time_range.sleepPRE(2,k) & event_time > time_range.sleepPRE(1,k);
                        %                             behavioural_state = -1;
                        %                             log_odd.behavioural_state(c)= behavioural_state;  % sleep pre
                        %                         end
                        %                     end
                        %
                        %                     if behavioural_state ~=-1
                        %                         for k = 1:size(time_range.awakePRE,2)
                        %                             if event_time < time_range.awakePRE(2,k) & event_time > time_range.awakePRE(1,k)
                        %                                 log_odd.behavioural_state(c)= 0;  % awake pre
                        %                             end
                        %                         end
                        %                     end
                        %
                        %                 elseif event_time>time_range.post(1) %(check just awake, sleep,seperated)
                        %                     for k = 1:size(time_range.sleepPOST,2)
                        %                         if event_time < time_range.sleepPOST(2,k) & event_time > time_range.sleepPOST(1,k)
                        %                             behavioural_state = 3;
                        %                             log_odd.behavioural_state(c)= behavioural_state;  % sleep POST
                        %                         end
                        %                     end
                        %
                        %                     if behavioural_state ~= 3
                        %                         for k = 1:size(time_range.awakePOST,2)
                        %                             if event_time < time_range.awakePOST(2,k) & event_time > time_range.awakePOST(1,k)
                        %                                 log_odd.behavioural_state(c)= 4;  % awake POST
                        %                             end
                        %                         end
                        %                     end
                        %
                        %                 elseif event_time>time_range.track(1).behaviour(1) & event_time<time_range.track(1).behaviour(2)
                        %                     log_odd.behavioural_state(c)=1;  % Run T1
                        %                 elseif event_time>time_range.track(2).behaviour(1) & event_time<time_range.track(2).behaviour(2)
                        %                     log_odd.behavioural_state(c)=2;  %Run T2
                        
                        
                    else
                        log_odd.behavioural_state(c)=NaN;  %boundary between behavioural states
                    end
                    
                    
                    % The 'ground truth' track ID label
                    log_odd.track(c)=track;
                    
                    % The session ID
                    log_odd.experiment(c)=f;
                    c=c+1;
                end
            end
        end

        % For multi-track events
        for i = 1:length(significant_replay_events.track(1).ref_index_MULTI)
            % save replay event index
            log_odd.index(c) = significant_replay_events.track(1).ref_index_MULTI(i);
            log_odd.multi_track(c) = 1;

            % The 'ground truth' track ID label (Find the track with better p val as the track id)
            event_index1 = find(significant_replay_events.track(1).ref_index == significant_replay_events.track(1).ref_index_MULTI(i));
            event_index2 = find(significant_replay_events.track(2).ref_index == significant_replay_events.track(1).ref_index_MULTI(i));

            if significant_replay_events.track(1).p_value(event_index1) < significant_replay_events.track(2).p_value(event_index2);
                track = 1;
            elseif significant_replay_events.track(2).p_value(event_index2) < significant_replay_events.track(1).p_value(event_index1);
                track = 2;
            else % if the p value is the same for both tracks,randomly assign a track id
                track = randperm(2);
                track = track(1);
            end

            if track == 1
                event_index = event_index1;
            elseif track == 2
                event_index = event_index2;
            end

            log_odd.track(c)=track;
            %             log_odd.pvalue(c)=significant_replay_events.track(track).p_value(event_index);
            log_odd.pvalue(c) = sig_event_info.p_value(track,significant_replay_events.track(track).index(event_index));
            log_odd.segment_id(c) = sig_event_info.segment_id(track,significant_replay_events.track(track).index(event_index));

            log_odd.neuron_count(c)= replay.neuron_count(log_odd.index(c));
            log_odd.spike_count(c)= replay.spike_count(log_odd.index(c));
 
%             log_odd.track1_pvalue(c)=log10(9e-4+ max(p_values.WHOLE(1,log_odd.index(c),:)));


            log_odd.best_score(c)=significant_replay_events.track(track).replay_score(event_index);
%             log_odd.track1_best_score(c)=replay_scores.WHOLE(1,log_odd.index(c));


            % Original 20ms
            log_odd.normal.probability_ratio(c) = probability_ratio_original{1}(track,log_odd.index(c),log_odd.segment_id(c));
            log_odd.normal.T1_T2_ratio(c) = probability_ratio_original{1}(1,log_odd.index(c),log_odd.segment_id(c));

            % True/False Ratio
            for nshuffle = 1:1000
                log_odd.normal.probability_ratio_shuffled{c}(nshuffle) = probability_ratio_original{2}{nshuffle}(track,log_odd.index(c),log_odd.segment_id(c));
            end

            % Calculate and save zscored log odd
            data = log(log_odd.normal.probability_ratio(c));
            shuffled_data = log(log_odd.normal.probability_ratio_shuffled{c});
            tempt = (data-mean(shuffled_data))/std(shuffled_data);
            log_odd.normal_zscored.original(1,c) = tempt;

            % T1/T2 ratio
            for nshuffle = 1:1000
                log_odd.normal.T1_T2_ratio_shuffled{c}(nshuffle) = probability_ratio_original{2}{nshuffle}(1,log_odd.index(c),log_odd.segment_id(c));
            end

            % Calculate and save zscored log odd
            data = log(log_odd.normal.T1_T2_ratio(c));
            shuffled_data = log(log_odd.normal.T1_T2_ratio_shuffled{c});
            tempt = (data-mean(shuffled_data))/std(shuffled_data);
            log_odd.normal_zscored.original(2,c) = tempt;


            % Original 20ms
            log_odd.common.probability_ratio(c) = probability_ratio_common_good{1}(track,log_odd.index(c),log_odd.segment_id(c));
            log_odd.common.T1_T2_ratio(c) = probability_ratio_common_good{1}(1,log_odd.index(c),log_odd.segment_id(c));

            % True/False Ratio
            for nshuffle = 1:1000
                log_odd.common.probability_ratio_shuffled{c}(nshuffle) = probability_ratio_common_good{2}{nshuffle}(track,log_odd.index(c),log_odd.segment_id(c));
            end

            % Calculate and save zscored log odd
            data = log(log_odd.common.probability_ratio(c));
            shuffled_data = log(log_odd.common.probability_ratio_shuffled{c});
            tempt = (data-mean(shuffled_data))/std(shuffled_data);
            log_odd.common_zscored.original(1,c) = tempt;

            % T1/T2 ratio
            for nshuffle = 1:1000
                log_odd.common.T1_T2_ratio_shuffled{c}(nshuffle) = probability_ratio_common_good{2}{nshuffle}(1,log_odd.index(c),log_odd.segment_id(c));
            end

            % Calculate and save zscored log odd
            data = log(log_odd.common.T1_T2_ratio(c));
            shuffled_data = log(log_odd.common.T1_T2_ratio_shuffled{c});
            tempt = (data-mean(shuffled_data))/std(shuffled_data);
            log_odd.common_zscored.original(2,c) = tempt;

            %                     log_odd.active_cells(c)= length(decoded_track_replay(1).replay_active_cells{log_odd.index(c)});
            log_odd.event_duration(c) = significant_replay_events.track(track).event_duration(event_index);
            log_odd.ripple_peak(c) = replay.ripple_peak(log_odd.index(c));

            if isempty(jump_distance{f}{track}{log_odd.index(c)})
                log_odd.max_jump(c) = nan;
            else
                log_odd.max_jump(c) = jump_distance{f}{track}{log_odd.index(c)};
            end

            event_time = (replay.onset(log_odd.index(c)) + replay.offset(log_odd.index(c)))/2; % Median event Time

            if event_time <= time_range.pre(2) %If PRE

                log_odd.behavioural_state(c)= -1;

            elseif event_time>=time_range.post(1) %If POST

                log_odd.behavioural_state(c)= 0;


            elseif event_time>=time_range.track(1).behaviour(1) & event_time<=time_range.track(1).behaviour(2)
                log_odd.behavioural_state(c)=1;  % If Run T1
            elseif event_time>=time_range.track(2).behaviour(1) & event_time<=time_range.track(2).behaviour(2)
                log_odd.behavioural_state(c)=2;  % If Run T2

            else
                log_odd.behavioural_state(c)=NaN;  %boundary between behavioural states
            end

            % The session ID
            log_odd.experiment(c)=f;
            c=c+1;
        end

    end
    
    
    if strcmp(option,'global remapped')
        cd .\global_remapped_shuffles
        DIR = dir('shuffle_*');
        DataPath = natsortfiles({DIR.name})'
        cd ..
        
        % for nfolders = 1:5
        for nfolders = 1:length(DataPath)
            destination = ['global_remapped_shuffles\shuffle_' num2str(nfolders)]
            copyfile('extracted_position.mat',destination)
            copyfile('extracted_sleep_state.mat',destination)
            cd(destination)

            load scored_replay
            load scored_replay_segments
%             method = 'wcorr';
            
            % gather p value and replay score information
            [p_values, replay_scores] = extract_score_and_pvalue_ground_truth(scored_replay, scored_replay1, scored_replay2, method, []);

            [significant_replay_events sig_event_info] = number_of_significant_replays_ground_truth(0.05,ripple_zscore_threshold,method,shuffles,[]);


%             if strcmp(method,'wcorr')
%                 save significant_replay_events_wcorr significant_replay_events
%             elseif strcmp(method,'spearman')
%                 save significant_replay_events_spearman significant_replay_events
%             elseif strcmp(method,'linear')
%                 save significant_replay_events_linear significant_replay_events
%             elseif strcmp(method,'path')
%                 save significant_replay_events_path significant_replay_events                
%             end
            % gather time information for sorting the event as PRE, POST or RUN
            [~,time_range]=sort_replay_events([],'wcorr');
            
            cd ..
            cd ..
            
            for track=1:2
                for i=1:length(significant_replay_events.track(track).p_value)
                    if ~ismember(significant_replay_events.track(track).ref_index(i),significant_replay_events.track(1).ref_index_MULTI) % if not multi-track
                        % save replay event index
                        log_odd.index(c) = significant_replay_events.track(track).ref_index(i);
                        log_odd.multi_track(c) = 0;

                        log_odd.pvalue(c) = sig_event_info.p_value(track,significant_replay_events.track(track).index(i));
                        log_odd.segment_id(c) = sig_event_info.segment_id(track,significant_replay_events.track(track).index(i));
                        log_odd.neuron_count(c)= replay.neuron_count(log_odd.index(c));
                        log_odd.spike_count(c)= replay.spike_count(log_odd.index(c));


                        % Original 20ms all good cells

                        log_odd.normal.global_remapped_probability_ratio(c) = global_remapped_original_probability_ratio{nfolders}{1}(track,log_odd.index(c),log_odd.segment_id(c));
                        log_odd.normal.global_remapped_T1_T2_ratio(c) = global_remapped_original_probability_ratio{nfolders}{1}(1,log_odd.index(c),log_odd.segment_id(c));

                        % True/False Ratio
                        for nshuffle = 1:1000
                            log_odd.normal.global_remapped_probability_ratio_shuffled{c}(nshuffle) = global_remapped_original_probability_ratio{nfolders}{2}{nshuffle}(track,log_odd.index(c),log_odd.segment_id(c));
                        end
                        
                        % Calculate and save zscored log odd
                        data = log(log_odd.normal.global_remapped_probability_ratio(c));
                        shuffled_data = log(log_odd.normal.global_remapped_probability_ratio_shuffled{c});
                        tempt = (data-mean(shuffled_data))/std(shuffled_data);
                        log_odd.normal_zscored.global_remapped_original(1,c) = tempt;
                        
                        % T1/T2 ratio
                        for nshuffle = 1:1000
                            log_odd.normal.global_remapped_T1_T2_ratio_shuffled{c}(nshuffle) = global_remapped_original_probability_ratio{nfolders}{2}{nshuffle}(1,log_odd.index(c),log_odd.segment_id(c));
                        end
                        
                        % Calculate and save zscored log odd
                        data = log(log_odd.normal.global_remapped_T1_T2_ratio(c));
                        shuffled_data = log(log_odd.normal.global_remapped_T1_T2_ratio_shuffled{c});
                        tempt = (data-mean(shuffled_data))/std(shuffled_data);
                        log_odd.normal_zscored.global_remapped_original(2,c) = tempt;
                        
                        
                        % global remapped common good
                        log_odd.common.global_remapped_probability_ratio(c) = global_remapped_common_good_probability_ratio{nfolders}{1}(track,log_odd.index(c),log_odd.segment_id(c));
                        log_odd.common.global_remapped_T1_T2_ratio(c) = global_remapped_common_good_probability_ratio{nfolders}{1}(1,log_odd.index(c),log_odd.segment_id(c));
                        
                        % True/False Ratio
                        for nshuffle = 1:1000
                            log_odd.common.global_remapped_probability_ratio_shuffled{c}(nshuffle) = global_remapped_common_good_probability_ratio{nfolders}{2}{nshuffle}(track,log_odd.index(c),log_odd.segment_id(c));
                        end
                        
                        % Calculate and save zscored log odd
                        data = log(log_odd.common.global_remapped_probability_ratio(c));
                        shuffled_data = log(log_odd.common.global_remapped_probability_ratio_shuffled{c});
                        tempt = (data-mean(shuffled_data))/std(shuffled_data);
                        log_odd.common_zscored.global_remapped_original(1,c) = tempt;
                        
                        % T1/T2 ratio
                        for nshuffle = 1:1000
                            log_odd.common.global_remapped_T1_T2_ratio_shuffled{c}(nshuffle) = global_remapped_common_good_probability_ratio{nfolders}{2}{nshuffle}(1,log_odd.index(c),log_odd.segment_id(c));
                        end
                        
                        % Calculate and save zscored log odd
                        data = log(log_odd.common.global_remapped_T1_T2_ratio(c));
                        shuffled_data = log(log_odd.common.global_remapped_T1_T2_ratio_shuffled{c});
                        tempt = (data-mean(shuffled_data))/std(shuffled_data);
                        log_odd.common_zscored.global_remapped_original(2,c) = tempt;


                        if isempty(jump_distance{nfolders}{f}{track}{log_odd.index(c)})
                            log_odd.max_jump(c) = nan;
                        else
                            log_odd.max_jump(c) = jump_distance{nfolders}{f}{track}{log_odd.index(c)};
                        end
                        
%                         log_odd.track1_pvalue(c)=log10(9e-4+ max(p_values.WHOLE(1,log_odd.index(c),:)));
                        log_odd.best_score(c)=significant_replay_events.track(track).replay_score(i);
%                         log_odd.track1_best_score(c)=replay_scores.WHOLE(1,log_odd.index(c));
                        %                 log_odd.active_cells(c)= length(decoded_track_replay(1).replay_active_cells{log_odd.index(c)});
                        log_odd.event_duration(c) = significant_replay_events.track(track).event_duration(i);
                        
                        event_time = (replay.onset(log_odd.index(c)) + replay.offset(log_odd.index(c)))/2; % Median event Time
                        
                        if event_time <= time_range.pre(2) %If PRE
                            
                            log_odd.behavioural_state(c)= -1;
                            
                        elseif event_time>=time_range.post(1) %If POST
                            
                            log_odd.behavioural_state(c)= 0;
                            
                            
                        elseif event_time>=time_range.track(1).behaviour(1) & event_time<=time_range.track(1).behaviour(2)
                            log_odd.behavioural_state(c)=1;  % If Run T1
                        elseif event_time>time_range.track(2).behaviour(1) && event_time<time_range.track(2).behaviour(2)
                            log_odd.behavioural_state(c)=2;  % If Run T2
                            
                            
                        else
                            log_odd.behavioural_state(c)=NaN;  %boundary between behavioural states
                        end
                        
                        
                        % The 'ground truth' track ID label
                        log_odd.track(c)=track;
                        
                        % The session ID
                        log_odd.experiment(c)=f;
                        c=c+1;

                    end
                end
            end


            % For multi-track events
            for i = 1:length(significant_replay_events.track(1).ref_index_MULTI)
                % save replay event index
                log_odd.index(c) = significant_replay_events.track(1).ref_index_MULTI(i);
                log_odd.multi_track(c) = 1;

                % The 'ground truth' track ID label (Find the track with better p val as the track id)
                event_index1 = find(significant_replay_events.track(1).ref_index == significant_replay_events.track(1).ref_index_MULTI(i));
                event_index2 = find(significant_replay_events.track(2).ref_index == significant_replay_events.track(1).ref_index_MULTI(i));

                if significant_replay_events.track(1).p_value(event_index1) < significant_replay_events.track(2).p_value(event_index2);
                    track = 1;
                elseif significant_replay_events.track(2).p_value(event_index2) < significant_replay_events.track(1).p_value(event_index1);
                    track = 2;
                else % if the p value is the same for both tracks,
                    track = randperm(2);
                    track = track(1);
                end

                if track == 1
                    event_index = event_index1;
                elseif track == 2
                    event_index = event_index2;
                end

                log_odd.track(c)=track;
                %                 log_odd.pvalue(c)=significant_replay_events.track(track).p_value(event_index);
                log_odd.pvalue(c) = sig_event_info.p_value(track,significant_replay_events.track(track).index(event_index));
                log_odd.segment_id(c) = sig_event_info.segment_id(track,significant_replay_events.track(track).index(event_index));

                log_odd.neuron_count(c)= replay.neuron_count(log_odd.index(c));
                log_odd.spike_count(c)= replay.spike_count(log_odd.index(c));

                %                 log_odd.track1_pvalue(c)=log10(9e-4+ max(p_values.WHOLE(1,log_odd.index(c),:)));
                log_odd.best_score(c)=significant_replay_events.track(track).replay_score(event_index);
                %                 log_odd.track1_best_score(c)=replay_scores.WHOLE(1,log_odd.index(c));

                % Global remapped original all good cells

                log_odd.normal.global_remapped_probability_ratio(c) = global_remapped_original_probability_ratio{nfolders}{1}(track,log_odd.index(c),log_odd.segment_id(c));
                log_odd.normal.global_remapped_T1_T2_ratio(c) = global_remapped_original_probability_ratio{nfolders}{1}(1,log_odd.index(c),log_odd.segment_id(c));

                % True/False Ratio
                for nshuffle = 1:1000
                    log_odd.normal.global_remapped_probability_ratio_shuffled{c}(nshuffle) = global_remapped_original_probability_ratio{nfolders}{2}{nshuffle}(track,log_odd.index(c),log_odd.segment_id(c));
                end

                % Calculate and save zscored log odd
                data = log(log_odd.normal.global_remapped_probability_ratio(c));
                shuffled_data = log(log_odd.normal.global_remapped_probability_ratio_shuffled{c});
                tempt = (data-mean(shuffled_data))/std(shuffled_data);
                log_odd.normal_zscored.global_remapped_original(1,c) = tempt;

                % T1/T2 ratio
                for nshuffle = 1:1000
                    log_odd.normal.global_remapped_T1_T2_ratio_shuffled{c}(nshuffle) = global_remapped_original_probability_ratio{nfolders}{2}{nshuffle}(1,log_odd.index(c),log_odd.segment_id(c));
                end

                % Calculate and save zscored log odd
                data = log(log_odd.normal.global_remapped_T1_T2_ratio(c));
                shuffled_data = log(log_odd.normal.global_remapped_T1_T2_ratio_shuffled{c});
                tempt = (data-mean(shuffled_data))/std(shuffled_data);
                log_odd.normal_zscored.global_remapped_original(2,c) = tempt;

                % global remapped common good
                log_odd.common.global_remapped_probability_ratio(c) = global_remapped_common_good_probability_ratio{nfolders}{1}(track,log_odd.index(c),log_odd.segment_id(c));
                log_odd.common.global_remapped_T1_T2_ratio(c) = global_remapped_common_good_probability_ratio{nfolders}{1}(1,log_odd.index(c),log_odd.segment_id(c));

                % True/False Ratio
                for nshuffle = 1:1000
                    log_odd.common.global_remapped_probability_ratio_shuffled{c}(nshuffle) = global_remapped_common_good_probability_ratio{nfolders}{2}{nshuffle}(track,log_odd.index(c),log_odd.segment_id(c));
                end

                % Calculate and save zscored log odd
                data = log(log_odd.common.global_remapped_probability_ratio(c));
                shuffled_data = log(log_odd.common.global_remapped_probability_ratio_shuffled{c});
                tempt = (data-mean(shuffled_data))/std(shuffled_data);
                log_odd.common_zscored.global_remapped_original(1,c) = tempt;

                % T1/T2 ratio
                for nshuffle = 1:1000
                    log_odd.common.global_remapped_T1_T2_ratio_shuffled{c}(nshuffle) = global_remapped_common_good_probability_ratio{nfolders}{2}{nshuffle}(1,log_odd.index(c),log_odd.segment_id(c));
                end

                % Calculate and save zscored log odd
                data = log(log_odd.common.global_remapped_T1_T2_ratio(c));
                shuffled_data = log(log_odd.common.global_remapped_T1_T2_ratio_shuffled{c});
                tempt = (data-mean(shuffled_data))/std(shuffled_data);
                log_odd.common_zscored.global_remapped_original(2,c) = tempt;




                %                     log_odd.active_cells(c)= length(decoded_track_replay(1).replay_active_cells{log_odd.index(c)});
                log_odd.event_duration(c) = significant_replay_events.track(track).event_duration(event_index);
                log_odd.ripple_peak(c) = replay.ripple_peak(log_odd.index(c));

                if isempty(jump_distance{nfolders}{f}{track}{log_odd.index(c)})
                    log_odd.max_jump(c) = nan;
                else
                    log_odd.max_jump(c) = jump_distance{nfolders}{f}{track}{log_odd.index(c)};
                end

                event_time = (replay.onset(log_odd.index(c)) + replay.offset(log_odd.index(c)))/2; % Median event Time

                if event_time <= time_range.pre(2) %If PRE
                    log_odd.behavioural_state(c)= -1;

                elseif event_time>=time_range.post(1) %If POST
                    log_odd.behavioural_state(c)= 0;

                elseif event_time>=time_range.track(1).behaviour(1) & event_time<=time_range.track(1).behaviour(2)
                    log_odd.behavioural_state(c)=1;  % If Run T1

                elseif event_time>=time_range.track(2).behaviour(1) & event_time<=time_range.track(2).behaviour(2)
                    log_odd.behavioural_state(c)=2;  % If Run T2

                else
                    log_odd.behavioural_state(c)=NaN;  %boundary between behavioural states
                end

                % The session ID
                log_odd.experiment(c)=f;
                c=c+1;
            end

        end
    end
    
    cd(current_directory)
end

end


