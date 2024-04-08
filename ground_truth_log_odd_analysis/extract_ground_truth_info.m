function [log_odd] = extract_ground_truth_info(folders,option,method,shuffles,ripple_zscore_threshold,p_val_threshold)
% Input:
% Folders: Each Folder is data for a session
% method = 'linear' or 'wcorr' or 'spearman' or
% p_val_threshold = If you want to extract events at specific p value, 
% set the p value threshold you want (e.g. 0.05) However, if you want to 
% extract and analyse all events, please set at 2 such that all events will be
% selected. (set at 2 such that HALF events with p value of 1 can be included for analysis) 

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
elseif strcmp(option,'global_remapped')
    load jump_distance_global_remapped
elseif strcmp(option,'spike_train_shifted')
    load jump_distance_spike_train_shifted
elseif strcmp(option,'place_field_shifted')
    load jump_distance_place_field_shifted
elseif strcmp(option,'cross_experiment_shuffled')
    load jump_distance_cross_experiment_shuffled
end
cd ..

for f = 1:length(folders)

    cd(folders{f})
    parameters = list_of_parameters;
    place_fields_BAYESIAN = [];
    %     load('significant_replay_events_wcorr.mat');
%     load('extracted_place_fields_BAYESIAN.mat');
    load('extracted_position.mat')
    load('extracted_clusters.mat');
    load('extracted_replay_events.mat');

    % Stable common good cells (based on subsets_of_cells)
    %     place_cell_index = subset_of_cells.cell_IDs{8}...
    %         (intersect(find(subset_of_cells.cell_IDs{8} >= 1000*f),find(subset_of_cells.cell_IDs{8} <= 1000*(f+1))))- f*1000;

    if strcmp(option,'original') || strcmp(option,'common')
        load probability_ratio_original
        load probability_ratio_common_good
    else % If randomised data, initialise the variables based on the shuffle name
        log_odd.normal.([option,'_T1_T2_ratio']) = [];
        log_odd.normal_zscored.([option,'_original']) = [];
        log_odd.common.([option,'_T1_T2_ratio']) = [];
        log_odd.common_zscored.([option,'_original']) = [];
    end

    if strcmp(option,'global_remapped')
        data_folders = 'global_remapped';
        load global_remapped_original_probability_ratio
        load global_remapped_common_good_probability_ratio

        probability_ratio_original = global_remapped_original_probability_ratio;
        probability_ratio_common_good = global_remapped_common_good_probability_ratio;
        clear global_remapped_original_probability_ratio global_remapped_common_good_probability_ratio
    elseif strcmp(option,'place_field_shifted')
        data_folders = 'place_field';
        load place_field_shifted_original_probability_ratio
        load place_field_shifted_common_good_probability_ratio

        probability_ratio_original = place_field_shifted_original_probability_ratio;
        probability_ratio_common_good = place_field_shifted_common_good_probability_ratio;
        clear place_field_shifted_original_probability_ratio place_field_shifted_common_good_probability_ratio
    elseif strcmp(option,'spike_train_shifted')
        data_folders = 'spike_train';
        load('spike_train_shifted_original_probability_ratio.mat')
        load spike_train_shifted_common_good_probability_ratio

        probability_ratio_original = spike_train_shifted_original_probability_ratio;
        probability_ratio_common_good = spike_train_shifted_common_good_probability_ratio;
        clear spike_train_shifted_original_probability_ratio spike_train_shifted_common_good_probability_ratio
    elseif strcmp(option,'cross_experiment_shuffled')
        data_folders = 'cross_experiment';
        load cross_experiment_shuffled_original_probability_ratio
        load cross_experiment_shuffled_common_good_probability_ratio

        probability_ratio_original = cross_experiment_shuffled_original_probability_ratio;
        probability_ratio_common_good = cross_experiment_shuffled_common_good_probability_ratio;
        clear cross_experiment_shuffled_original_probability_ratio cross_experiment_shuffled_common_good_probability_ratio
        
    end



    %% Ground Truth Information (Bayesian bias, p value, best score, probability ratio)


    if strcmp(option,'original') || strcmp(option,'common')

        load scored_replay
        load scored_replay_segments

        % gather p value and replay score information
        [p_values, replay_scores] = extract_score_and_pvalue_ground_truth(scored_replay, scored_replay1, scored_replay2, method, []);

        [significant_replay_events sig_event_info] = number_of_significant_replays_ground_truth(p_val_threshold,ripple_zscore_threshold,method,shuffles,[]);

        % gather time information for sorting the event as PRE, POST or RUN
        [~,time_range]=sort_replay_events([],'wcorr');

        for i=1:length(significant_replay_events.track(1).p_value)
            for track=1:2
                % save replay event index
                log_odd.index(c) = significant_replay_events.track(track).ref_index(i);

                log_odd.pvalue(track,c) = sig_event_info.p_value(track,significant_replay_events.track(track).index(i));
                log_odd.segment_id(track,c) = sig_event_info.segment_id(track,significant_replay_events.track(track).index(i));

                % Original 20ms
                %                 log_odd.normal.probability_ratio(track,c) = probability_ratio_original{1}(track,log_odd.index(c),log_odd.segment_id(c));
                log_odd.normal.T1_T2_ratio(track,c) = probability_ratio_original{1}(1,log_odd.index(c),log_odd.segment_id(track,c));

                %                 % True/False Ratio
                %                 for nshuffle = 1:1000
                %                     log_odd.normal.probability_ratio_shuffled{c}(nshuffle) = probability_ratio_original{2}{nshuffle}(track,log_odd.index(c),log_odd.segment_id(c));
                %                 end
                %
                %                 % Calculate and save zscored log odd
                %                 data = log(log_odd.normal.probability_ratio(track,c));
                %                 shuffled_data = log(log_odd.normal.probability_ratio_shuffled{track}{c});
                %                 tempt = (data-mean(shuffled_data))/std(shuffled_data);
                %                 log_odd.normal_zscored.original(1,c) = tempt;
                T1_T2_ratio_shuffled = [];
                % T1/T2 ratio
                for nshuffle = 1:1000
                    T1_T2_ratio_shuffled(nshuffle) = probability_ratio_original{2}{nshuffle}(1,log_odd.index(c),log_odd.segment_id(track,c));
                end

                % Calculate and save zscored log odd
                data = log(log_odd.normal.T1_T2_ratio(track,c));
                shuffled_data = log(T1_T2_ratio_shuffled);
                tempt = (data-mean(shuffled_data))/std(shuffled_data);
                log_odd.normal_zscored.original(track,c) = tempt;


                % Common 20ms
                %                 log_odd.common.probability_ratio(c) = probability_ratio_common_good{1}(track,log_odd.index(c),log_odd.segment_id(c));
                log_odd.common.T1_T2_ratio(track,c) = probability_ratio_common_good{1}(1,log_odd.index(c),log_odd.segment_id(track,c));

                %                 % True/False Ratio
                %                 for nshuffle = 1:1000
                %                     log_odd.common.probability_ratio_shuffled{c}(nshuffle) = probability_ratio_common_good{2}{nshuffle}(track,log_odd.index(c),log_odd.segment_id(c));
                %                 end
                %
                %                 % Calculate and save zscored log odd
                %                 data = log(log_odd.common.probability_ratio(c));
                %                 shuffled_data = log(log_odd.common.probability_ratio_shuffled{c});
                %                 tempt = (data-mean(shuffled_data))/std(shuffled_data);
                %                 log_odd.common_zscored.original(1,c) = tempt;
                T1_T2_ratio_shuffled = [];
                % T1/T2 ratio
                for nshuffle = 1:1000
                    T1_T2_ratio_shuffled(nshuffle) = probability_ratio_common_good{2}{nshuffle}(1,log_odd.index(c),log_odd.segment_id(track,c));
                end

                % Calculate and save zscored log odd
                data = log(log_odd.common.T1_T2_ratio(track,c));
                shuffled_data = log(T1_T2_ratio_shuffled);
                tempt = (data-mean(shuffled_data))/std(shuffled_data);
                log_odd.common_zscored.original(track,c) = tempt;


                %                     log_odd.track1_pvalue(c)=log10(9e-4+  max(p_values.WHOLE(1,log_odd.index(c),:)));
                %                     log_odd.track1_pvalue(c)=log10(sig_event_info.p_value(1,log_odd.index(c)));
                %                     log_odd.pvalue(c)=log10(9e-4+ significant_replay_events.track(track).p_value(i));


                log_odd.best_score(track,c)=significant_replay_events.track(track).replay_score(i);
                %                     log_odd.track1_best_score(c)=replay_scores.WHOLE(1,log_odd.index(c));

                log_odd.neuron_count(c)= replay.neuron_count(log_odd.index(c));
                log_odd.spike_count(c)= replay.spike_count(log_odd.index(c));

                log_odd.event_duration(c) = significant_replay_events.track(track).event_duration(i);
                log_odd.ripple_peak(c) = replay.ripple_peak(log_odd.index(c));

                if isempty(jump_distance{f}{track}{log_odd.index(c)})
                    log_odd.max_jump(track,c) = nan;
                else
                    log_odd.max_jump(track,c) = jump_distance{f}{track}{log_odd.index(c)};
                end


                event_time = (replay.onset(log_odd.index(c)) + replay.offset(log_odd.index(c)))/2; % Median event Time

                %                     event_time = significant_replay_events.track(track).event_times(i);

                if event_time >= time_range.pre(1) & event_time <= time_range.pre(2) %If PRE

                    log_odd.behavioural_state(c)= -1;

                elseif event_time >= time_range.post(1) & event_time <= time_range.post(2) %If POST

                    log_odd.behavioural_state(c)= 0;


                elseif event_time>=time_range.track(1).behaviour(1) & event_time<=time_range.track(1).behaviour(2)
                    log_odd.behavioural_state(c)=1;  % If Run T1
                elseif event_time>=time_range.track(2).behaviour(1) & event_time<=time_range.track(2).behaviour(2)
                    log_odd.behavioural_state(c)=2;  % If Run T2

                else
                    log_odd.behavioural_state(c)=NaN;  %boundary between behavioural states
                end


                % The 'ground truth' track ID label
                log_odd.track(track,c)=track;

                % The session ID
                log_odd.experiment(c)=f;
            end
            c=c+1;
        end


    else % Randomised dataset
        
        cd([data_folders,'_shuffles'])
        DIR = dir('shuffle_*');
        DataPath = natsortfiles({DIR.name})'
        cd ..

        % gather time information for sorting the event as PRE, POST or RUN
        [~,time_range]=sort_replay_events([],'wcorr');

        % for nfolders = 1:5
        for nfolders = 1:length(DataPath)
            destination = [data_folders,'_shuffles\shuffle_' num2str(nfolders)]

            if isfile('extracted_sleep_state.mat')
                copyfile('extracted_replay_events.mat',destination)
                copyfile('extracted_position.mat',destination)
                copyfile('extracted_sleep_state.mat',destination)
            end

            cd(destination)
            load scored_replay
            load scored_replay_segments
            %             method = 'wcorr';

            % gather p value and replay score information
            [p_values, replay_scores] = extract_score_and_pvalue_ground_truth(scored_replay, scored_replay1, scored_replay2, method, []);

            [significant_replay_events sig_event_info] = number_of_significant_replays_ground_truth(p_val_threshold,ripple_zscore_threshold,method,shuffles,[]);

            cd ..
            cd ..

            for i=1:length(significant_replay_events.track(1).p_value)
                for track=1:2
                    % save replay event index
                    log_odd.index(c) = significant_replay_events.track(track).ref_index(i);

                    log_odd.pvalue(track,c) = sig_event_info.p_value(track,significant_replay_events.track(track).index(i));
                    log_odd.segment_id(track,c) = sig_event_info.segment_id(track,significant_replay_events.track(track).index(i));

                    % Original 20ms
                    %                 log_odd.normal.probability_ratio(track,c) = probability_ratio_original{1}(track,log_odd.index(c),log_odd.segment_id(c));
                    log_odd.normal.([option,'_T1_T2_ratio'])(track,c) = probability_ratio_original{nfolders}{1}(1,log_odd.index(c),log_odd.segment_id(track,c));

                    %                 % True/False Ratio
                    %                 for nshuffle = 1:1000
                    %                     log_odd.normal.probability_ratio_shuffled{c}(nshuffle) = probability_ratio_original{2}{nshuffle}(track,log_odd.index(c),log_odd.segment_id(c));
                    %                 end
                    %
                    %                 % Calculate and save zscored log odd
                    %                 data = log(log_odd.normal.probability_ratio(track,c));
                    %                 shuffled_data = log(log_odd.normal.probability_ratio_shuffled{track}{c});
                    %                 tempt = (data-mean(shuffled_data))/std(shuffled_data);
                    %                 log_odd.normal_zscored.original(1,c) = tempt;
                    T1_T2_ratio_shuffled = [];
                    % T1/T2 ratio
                    for nshuffle = 1:1000
                        T1_T2_ratio_shuffled(nshuffle) = probability_ratio_original{nfolders}{2}{nshuffle}(1,log_odd.index(c),log_odd.segment_id(track,c));
                    end

                    % Calculate and save zscored log odd
                    data = log(log_odd.normal.([option,'_T1_T2_ratio'])(track,c));
                    shuffled_data = log(T1_T2_ratio_shuffled);
                    tempt = (data-mean(shuffled_data))/std(shuffled_data);
                    log_odd.normal_zscored.([option,'_original'])(track,c) = tempt;


                    % Common 20ms
                    %                 log_odd.common.probability_ratio(c) = probability_ratio_common_good{1}(track,log_odd.index(c),log_odd.segment_id(c));
                    log_odd.common.([option,'_T1_T2_ratio'])(track,c) = probability_ratio_common_good{nfolders}{1}(1,log_odd.index(c),log_odd.segment_id(track,c));

                    %                 % True/False Ratio
                    %                 for nshuffle = 1:1000
                    %                     log_odd.common.probability_ratio_shuffled{c}(nshuffle) = probability_ratio_common_good{2}{nshuffle}(track,log_odd.index(c),log_odd.segment_id(c));
                    %                 end
                    %
                    %                 % Calculate and save zscored log odd
                    %                 data = log(log_odd.common.probability_ratio(c));
                    %                 shuffled_data = log(log_odd.common.probability_ratio_shuffled{c});
                    %                 tempt = (data-mean(shuffled_data))/std(shuffled_data);
                    %                 log_odd.common_zscored.original(1,c) = tempt;
                    T1_T2_ratio_shuffled = [];
                    % T1/T2 ratio
                    for nshuffle = 1:1000
                        T1_T2_ratio_shuffled(nshuffle) = probability_ratio_common_good{nfolders}{2}{nshuffle}(1,log_odd.index(c),log_odd.segment_id(track,c));
                    end

                    % Calculate and save zscored log odd
                    data = log(log_odd.common.([option,'_T1_T2_ratio'])(track,c));
                    shuffled_data = log(T1_T2_ratio_shuffled);
                    tempt = (data-mean(shuffled_data))/std(shuffled_data);
                    log_odd.common_zscored.([option,'_original'])(track,c) = tempt;


                    %                     log_odd.track1_pvalue(c)=log10(9e-4+  max(p_values.WHOLE(1,log_odd.index(c),:)));
                    %                     log_odd.track1_pvalue(c)=log10(sig_event_info.p_value(1,log_odd.index(c)));
                    %                     log_odd.pvalue(c)=log10(9e-4+ significant_replay_events.track(track).p_value(i));


                    log_odd.best_score(track,c)=significant_replay_events.track(track).replay_score(i);
                    %                     log_odd.track1_best_score(c)=replay_scores.WHOLE(1,log_odd.index(c));

                    log_odd.neuron_count(c)= replay.neuron_count(log_odd.index(c));
                    log_odd.spike_count(c)= replay.spike_count(log_odd.index(c));

                    log_odd.event_duration(c) = significant_replay_events.track(track).event_duration(i);
                    log_odd.ripple_peak(c) = replay.ripple_peak(log_odd.index(c));

                    if isempty(jump_distance{nfolders}{f}{track}{log_odd.index(c)})
                        log_odd.max_jump(track,c) = nan;
                    else
                        log_odd.max_jump(track,c) = jump_distance{nfolders}{f}{track}{log_odd.index(c)};
                    end


                    event_time = (replay.onset(log_odd.index(c)) + replay.offset(log_odd.index(c)))/2; % Median event Time

                    %                     event_time = significant_replay_events.track(track).event_times(i);

                    if event_time >= time_range.pre(1) & event_time <= time_range.pre(2) %If PRE

                        log_odd.behavioural_state(c)= -1;

                    elseif event_time >= time_range.post(1) & event_time <= time_range.post(2) %If POST

                        log_odd.behavioural_state(c)= 0;


                    elseif event_time>=time_range.track(1).behaviour(1) & event_time<=time_range.track(1).behaviour(2)
                        log_odd.behavioural_state(c)=1;  % If Run T1
                    elseif event_time>=time_range.track(2).behaviour(1) & event_time<=time_range.track(2).behaviour(2)
                        log_odd.behavioural_state(c)=2;  % If Run T2

                    else
                        log_odd.behavioural_state(c)=NaN;  %boundary between behavioural states
                    end


                    % The 'ground truth' track ID label
                    log_odd.track(track,c)=track;

                    % The session ID
                    log_odd.experiment(c)=f;

                    % The shuffle ID
                    log_odd.shuffle(c)=nfolders;
                end
                c=c+1;
            end

        end


    end

    cd(current_directory)

%     if strcmp(option,'original') || strcmp(option,'common')
%         %     log_odd.common = rmfield(log_odd.common,{'T1_T2_ratio_shuffled','probability_ratio_shuffled'});
%         %     log_odd.normal = rmfield(log_odd.normal,{'T1_T2_ratio_shuffled','probability_ratio_shuffled'});
%         log_odd.common = rmfield(log_odd.common,{'T1_T2_ratio_shuffled'});
%         log_odd.normal = rmfield(log_odd.normal,{'T1_T2_ratio_shuffled'});
%     else
%         log_odd.common = rmfield(log_odd.common,{'global_remapped_T1_T2_ratio_shuffled'});
%         log_odd.normal = rmfield(log_odd.normal,{'global_remapped_T1_T2_ratio_shuffled'});
%         %         log_odd.common = rmfield(log_odd.common,{'global_remapped_T1_T2_ratio_shuffled','global_remapped_probability_ratio_shuffled'});
%         %     log_odd.normal = rmfield(log_odd.normal,{'global_remapped_T1_T2_ratio_shuffled','global_remapped_probability_ratio_shuffled'});
%     end

clear probability_ratio_original probability_ratio_common_good

end

end


