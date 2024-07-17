function [] = plot_ground_truth_compare_multitrack(folders,option,method)

total_number{1}(1:3) = 0;
total_number{2}(1:3) = 0;

for nfolder = 1:length(folders) % Get total number of replay events
    for nshuffle = 1:2 % 1 is original and 2 is global remapped
        session_total_number{nshuffle}(nfolder,1:3) = 0;
    end
end

for nfolder = 1:length(folders) % Get total number of replay events
    for nshuffle = 1:2 % 1 is original and 2 is global remapped

        cd(folders{nfolder})

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

            elseif event_time >= time_range.track(1).behaviour(1) & event_time <= time_range.track(1).behaviour(2) % If RUN Track 1
                total_number{nshuffle}(2) = total_number{nshuffle}(2) + 1*number_of_global_remapped_shuffles;
                session_total_number{nshuffle}(nfolder,2) = session_total_number{nshuffle}(nfolder,2) + 1*number_of_global_remapped_shuffles;

            elseif event_time >= time_range.track(2).behaviour(1) & event_time <= time_range.track(2).behaviour(2) % If RUN Track 2
                total_number{nshuffle}(2) = total_number{nshuffle}(2) + 1*number_of_global_remapped_shuffles;
                session_total_number{nshuffle}(nfolder,2) = session_total_number{nshuffle}(nfolder,2) + 1*number_of_global_remapped_shuffles;

            end
        end
        cd ..
    end
end


index = [];
epoch_index = [];
for nmethod = 1:length(method)
    cd ground_truth_original
    if strcmp(method{nmethod},'wcorr 1 shuffle')
        load log_odd_wcorr_place_bin_circular_shift
        log_odd_compare{nmethod}{1} = log_odd;
        load log_odd_wcorr_place_bin_circular_shift_global_remapped
        log_odd_compare{nmethod}{2} = log_odd;

    elseif strcmp(method{nmethod},'wcorr 1 shuffle + jump distance')
        load log_odd_wcorr_place_bin_circular_shift
        log_odd_compare{nmethod}{1} = log_odd;
        load log_odd_wcorr_place_bin_circular_shift_global_remapped
        log_odd_compare{nmethod}{2} = log_odd;

    elseif strcmp(method{nmethod},'wcorr 2 shuffles')
        load log_odd_wcorr_PRE_place_POST_time
        log_odd_compare{nmethod}{1} = log_odd;
        load log_odd_wcorr_PRE_place_POST_time_global_remapped
        log_odd_compare{nmethod}{2} = log_odd;

    elseif strcmp(method{nmethod},'wcorr 3 shuffles')
        load log_odd_wcorr
        log_odd_compare{nmethod}{1} = log_odd;
        load log_odd_wcorr_global_remapped
        log_odd_compare{nmethod}{2} = log_odd;


    elseif strcmp(method{nmethod},'spearman median spike')
        load log_odd_spearman
        log_odd_compare{nmethod}{1} = log_odd;
        load log_odd_spearman_global_remapped
        log_odd_compare{nmethod}{2} = log_odd;
    elseif strcmp(method{nmethod},'spearman all spikes')
        load log_odd_spearman_all_spikes
        log_odd_compare{nmethod}{1} = log_odd;
        load log_odd_spearman_all_spikes_global_remapped
        log_odd_compare{nmethod}{2} = log_odd;
        
    elseif strcmp(method{nmethod},'linear 1 shuffle')
        load log_odd_linear_place_bin_circular_shift
        log_odd_compare{nmethod}{1} = log_odd;
        load log_odd_linear_place_bin_circular_shift_global_remapped
        log_odd_compare{nmethod}{2} = log_odd;
        
    elseif strcmp(method{nmethod},'linear 2 shuffles')
        load log_odd_linear_PRE_place_POST_time
        log_odd_compare{nmethod}{1} = log_odd;
        load log_odd_linear_PRE_place_POST_time_global_remapped
        log_odd_compare{nmethod}{2} = log_odd;
    end
    cd ..

    % Get index for PRE,  RUN, POST
    states = [-1 0 1 2]; % PRE, POST, RUN Track 1, RUN Track 2
    
    for nshuffle = 1:2
        for k=1:length(states) % sort track id and behavioral states (pooled from 5 sessions)
            % Find intersect of behavioural state and ripple peak threshold
            state_index = find(log_odd_compare{nmethod}{nshuffle}.behavioural_state==states(k));
            
            if k == 1 % PRE
                epoch_index{nmethod}{nshuffle}{1} = state_index;
            elseif k == 2 % POST
                epoch_index{nmethod}{nshuffle}{3} = state_index;
            elseif k == 3 % RUN Track 1
                epoch_index{nmethod}{nshuffle}{2} = state_index;
            elseif k == 4 % RUN Track 2
                epoch_index{nmethod}{nshuffle}{2} = [epoch_index{nmethod}{nshuffle}{2} state_index];
            end
        end
    end

    % Grab all the essential information including:
    % p value
    % segment with best p value for a given track (whole event, first half or second half)
    % maximum jump distance
    % replay score (e.g. weighted correlation score)

    if strcmp(option,'common')
        data{nmethod}{1} = log_odd_compare{nmethod}{1}.common_zscored.original;
        data{nmethod}{2} = log_odd_compare{nmethod}{2}.common_zscored.global_remapped_original;

        for nshuffle = 1:length(log_odd_compare{nmethod})
            log_pval{nmethod}{nshuffle} = log_odd_compare{nmethod}{nshuffle}.pvalue;
            segment_id{nmethod}{nshuffle} = log_odd_compare{nmethod}{nshuffle}.segment_id;
            max_jump{nmethod}{nshuffle} = log_odd_compare{nmethod}{nshuffle}.max_jump;
            replay_score{nmethod}{nshuffle} = log_odd_compare{nmethod}{nshuffle}.best_score;
        end

    elseif strcmp(option,'original')
        data{nmethod}{1} = log_odd_compare{nmethod}{1}.normal_zscored.original;
        data{nmethod}{2} = log_odd_compare{nmethod}{2}.normal_zscored.global_remapped_original;

        for nshuffle = 1:length(log_odd_compare{nmethod})
            log_pval{nmethod}{nshuffle} = log_odd_compare{nmethod}{nshuffle}.pvalue;
            segment_id{nmethod}{nshuffle} = log_odd_compare{nmethod}{nshuffle}.segment_id;
            max_jump{nmethod}{nshuffle} = log_odd_compare{nmethod}{nshuffle}.max_jump;
            replay_score{nmethod}{nshuffle} = log_odd_compare{nmethod}{nshuffle}.best_score;
        end
    end

    %     colour_line2 = {'k--','b--','r--','g--'};
    colour_symbol={'bo','ro','go','ko'};
    Behavioural_epoches = {'PRE','RUN','POST'};
end


p_val_threshold = round(10.^[-3:0.01:log10(0.2)],4);
p_val_threshold = flip(unique(p_val_threshold));
p_val_threshold(1) = 0.2;
low_threshold = find(ismembertol(p_val_threshold,[0.0501 0.02 0.01 0.005 0.002 0.001],eps) == 1);

cd ground_truth_original
if exist('log_odd_difference_optimisation.mat', 'file') ~= 2
    
    
    load log_odd_difference_multiple_shuffles
    % copy wcorr 1 shuffle
    log_odd_difference1{1} = log_odd_difference{1};
    log_odd_difference_CI1{1} = log_odd_difference_CI{1};
    percent_sig_events1{1} = percent_sig_events{1};
    percent_sig_events_CI1{1} = percent_sig_events_CI{1};
    percent_multi_events1{1} = percent_multi_events{1};
    percent_multi_events_CI1{1} = percent_multi_events_CI{1};
    percent_shuffle_events1{1} = percent_shuffle_events{1};
    percent_shuffle_events_CI1{1} = percent_shuffle_events_CI{1};
    replay_score_compare1{1} = replay_score_compare{1};
    replay_score_compare_CI1{1} = replay_score_compare_CI{1};

    % copy wcorr 2 shuffles
    log_odd_difference1{3} = log_odd_difference{4};
    log_odd_difference_CI1{3} = log_odd_difference_CI{4};
    percent_sig_events1{3} = percent_sig_events{4};
    percent_sig_events_CI1{3} = percent_sig_events_CI{4};
    percent_multi_events1{3} = percent_multi_events{4};
    percent_multi_events_CI1{3} = percent_multi_events_CI{4};
    percent_shuffle_events1{3} = percent_shuffle_events{4};
    percent_shuffle_events_CI1{3} = percent_shuffle_events_CI{4};
    replay_score_compare1{3} = replay_score_compare{4};
    replay_score_compare_CI1{3} = replay_score_compare_CI{4};

    % copy wcorr 3 shuffles
    log_odd_difference1{4} = log_odd_difference{7};
    log_odd_difference_CI1{4} = log_odd_difference_CI{7};
    percent_sig_events1{4} = percent_sig_events{7};
    percent_sig_events_CI1{4} = percent_sig_events_CI{7};
    percent_multi_events1{4} = percent_multi_events{7};
    percent_multi_events_CI1{4} = percent_multi_events_CI{7};
    percent_shuffle_events1{4} = percent_shuffle_events{7};
    percent_shuffle_events_CI1{4} = percent_shuffle_events_CI{7};
    replay_score_compare1{4} = replay_score_compare{7};
    replay_score_compare_CI1{4} = replay_score_compare_CI{7};

    load log_odd_difference_jump_distance
    % copy wcorr 1 shuffle + jump distance 0.4 (From 40 because it equals to p <= 0.02)
    for epoch = 1:3
        for nshuffle = 1:2
            log_odd_difference1{2}{nshuffle}{epoch} = log_odd_difference{3}{nshuffle}{epoch}(71:end,:);
            log_odd_difference_CI1{2}{nshuffle}{epoch} = log_odd_difference_CI{3}{nshuffle}{epoch}(71:end,:);
            percent_sig_events1{2}{nshuffle}{epoch} = percent_sig_events{3}{nshuffle}{epoch}(71:end,:);
            percent_sig_events_CI1{2}{nshuffle}{epoch} = percent_sig_events_CI{3}{nshuffle}{epoch}(71:end,:);
            percent_multi_events1{2}{nshuffle}{epoch} = percent_multi_events{3}{nshuffle}{epoch}(71:end,:);
            percent_multi_events_CI1{2}{nshuffle}{epoch} = percent_multi_events_CI{3}{nshuffle}{epoch}(71:end,:);
            replay_score_compare1{2}{nshuffle}{epoch} = replay_score_compare{3}{nshuffle}{epoch}(71:end,:)
            replay_score_compare_CI1{2}{nshuffle}{epoch} = replay_score_compare_CI{3}{nshuffle}{epoch}(71:end,:);
        end
        percent_shuffle_events1{2}{epoch} = percent_shuffle_events{3}{epoch}(71:end,:);
        percent_shuffle_events_CI1{2}{epoch} = percent_shuffle_events_CI{3}{epoch}(71:end,:);
    end
    
    % putative analysis of p value distribution
    for nmethod = 6
        for epoch = 1:3
            %                     % Resample canditate events with replacement
            %                     resampled_event = datasample(s,epoch_index{nmethod}{nshuffle}{epoch},length(epoch_index{nmethod}{nshuffle}{epoch}));
            %                     resampled_event_p_value = log_pval{nmethod}{nshuffle}(1,resampled_event);
            %
            figure(nmethod)
            sgtitle('Track 1')
            subplot(2,2,epoch)
            % Actual data
            event_p_value1 = log_pval{nmethod}{1}(1,epoch_index{nmethod}{1}{epoch});
           
           q_values = (1:length(event_p_value1))/length(event_p_value1)   * 0.05; % calculate critical values

            [sorted_p, sort_idx] = sort(event_p_value1); % sort p-values in ascending order

            % Reject all candidate events with p-values less than or equal to the q-value
            if ~isempty(find(sorted_p <= sorted_q, 1, 'last'))% find last significant event
                significant_idx = sort_idx(1:find(sorted_p <= sorted_q, 1, 'last'));
                significant_proportion = length(significant_idx)/length(event_p_value1);
            else
                % no significant events found
            end

             histogram(event_p_value1,0:0.02:1,'Normalization','probability','EdgeAlpha',0.1)
            hold on
            % cell-id shuffled data
            event_p_value2 = sort(log_pval{nmethod}{2}(1,epoch_index{nmethod}{2}{epoch}));
            histogram(event_p_value2,0:0.02:1,'Normalization','probability','EdgeAlpha',0.1)
            legend('Actual data','cell-id shuffled data')
            title(Behavioural_epoches{epoch})
            
            figure(length(log_pval)*2+nmethod)
            sgtitle('Track 1')
            subplot(2,2,epoch)
            [fdr,q,priori] = mafdr(event_p_value1,'Showplot',false)
%             plot(event_p_value1,q);
            plot(1:length(event_p_value1),event_p_value1);
            hold on
            [fdr,q,priori] = mafdr(event_p_value2,'Showplot',false)
%             plot(event_p_value2,q);
            s = RandStream('mrg32k3a','Seed',1); % Set random seed for resampling
            resampled_event = sort(datasample(s,event_p_value2,length(event_p_value1)));
            plot(1:length(resampled_event),resampled_event);

            figure(length(log_pval) + nmethod)
            sgtitle('Track 2')
            subplot(2,2,epoch)
            % Actual data
            event_p_value1 = sort(log_pval{nmethod}{1}(2,epoch_index{nmethod}{1}{epoch}));
             histogram(event_p_value1,0:0.02:1,'Normalization','probability','EdgeAlpha',0.1)
            hold on
            % cell-id shuffled data
            event_p_value2 = sort(log_pval{nmethod}{2}(2,epoch_index{nmethod}{2}{epoch}));
             histogram(event_p_value2,0:0.02:1,'Normalization','probability','EdgeAlpha',0.1)
            legend('Actual data','cell-id shuffled data')
            title(Behavioural_epoches{epoch})
            
            figure(length(log_pval)*3+nmethod)
            sgtitle('Track 2')
            subplot(2,2,epoch)
            [fdr,q,priori] = mafdr(event_p_value1,'Showplot',false)
%             plot(event_p_value1,q);
            plot(1:length(event_p_value1),event_p_value1);
            hold on
            [fdr,q,priori] = mafdr(event_p_value2,'Showplot',false)
            %             plot(event_p_value2,q);
            s = RandStream('mrg32k3a','Seed',1); % Set random seed for resampling
            resampled_event = sort(datasample(s,event_p_value2,length(event_p_value1)));
            plot(1:length(resampled_event),resampled_event);
     
        end
    end

    for nmethod = 5:length(method)
        for epoch = 1:3
            for nshuffle = 1:length(log_pval{nmethod})
                tic
                
                parfor threshold = 1:length(p_val_threshold)
                    for nboot = 1:1000 % Bootstrapping 1000 times
                        s = RandStream('mrg32k3a','Seed',nboot); % Set random seed for resampling

                        % Resample canditate events with replacement
                        resampled_event = datasample(s,epoch_index{nmethod}{nshuffle}{epoch},length(epoch_index{nmethod}{nshuffle}{epoch}));
                        
                        % Detecting events significant for track 1
                        % Find significant 'WHOLE' event
                        this_segment_event = resampled_event(find(segment_id{nmethod}{nshuffle}(1,resampled_event) == 1));
                        track_1_index = this_segment_event(find(log_pval{nmethod}{nshuffle}(1,this_segment_event) <= p_val_threshold(threshold)));
                        % Find significant 'half' event
                        this_segment_event = resampled_event(find(segment_id{nmethod}{nshuffle}(1,resampled_event) > 1));
                        track_1_index = [track_1_index this_segment_event(find(log_pval{nmethod}{nshuffle}(1,this_segment_event) <= p_val_threshold(threshold)/2))];% p value divide by two to account for multiple comparision

                        % Detecting events significant for track 2
                        % Find significant 'WHOLE' event
                        this_segment_event = resampled_event(find(segment_id{nmethod}{nshuffle}(2,resampled_event) == 1));
                        track_2_index = this_segment_event(find(log_pval{nmethod}{nshuffle}(2,this_segment_event) <= p_val_threshold(threshold)));
                        % Find significant 'half' event
                        this_segment_event = resampled_event(find(segment_id{nmethod}{nshuffle}(2,resampled_event) > 1));
                        track_2_index = [track_2_index this_segment_event(find(log_pval{nmethod}{nshuffle}(2,this_segment_event) <= p_val_threshold(threshold)/2))];% p value divide by two to account for multiple comparision

                        if strcmp(method{nmethod},'wcorr 1 shuffle + jump distance') % Apply jump distance threshold 0.4
                            track_1_index = track_1_index(find(max_jump{nmethod}{nshuffle}(1,track_1_index) <= 20*0.4));
                            track_2_index = track_2_index(find(max_jump{nmethod}{nshuffle}(2,track_2_index) <= 20*0.4));
                        end
                        
                        % Calculate mean replay score for all significant
                        % events
                        boot_replay_score(threshold,nboot) = mean([replay_score{nmethod}{nshuffle}(1,track_1_index) replay_score{nmethod}{nshuffle}(2,track_2_index)]);

                        % Calculate percentage of tracks significant for
                        % both tracks
                        multi_event_number = sum(ismember(track_1_index,track_2_index));
                        multi_event_percent = multi_event_number/(length(track_1_index)+length(track_2_index)-multi_event_number);
                        
                        % Calculate mean z-scored log odds difference for all significant
                        % events (Track 1 - Track 2)
                        boot_log_odd_difference(threshold,nboot) = mean(data{nmethod}{nshuffle}(1,track_1_index)) ...
                            - mean(data{nmethod}{nshuffle}(2,track_2_index));
                        

                        boot_percent_multi_events(threshold,nboot) = multi_event_percent;

                        % Significant event proportion (minus multitrack event number to avoid double counting)
                        boot_percent_sig_events(threshold,nboot) = (length(track_1_index)+length(track_2_index)-multi_event_number)/total_number{nshuffle}(epoch);

                        if nshuffle == 2
                            % cell id shuffled mean significant event proportion
                            boot_percent_shuffle_events(threshold,nboot) = ((length(track_1_index)+length(track_2_index))/2) / total_number{nshuffle}(epoch);
                        end
                    end
                end
                
                log_odd_difference1{nmethod}{nshuffle}{epoch} = boot_log_odd_difference;
                log_odd_difference_CI1{nmethod}{nshuffle}{epoch} = prctile(boot_log_odd_difference,[2.5 97.5],2);
                percent_sig_events1{nmethod}{nshuffle}{epoch} = boot_percent_sig_events;
                percent_sig_events_CI1{nmethod}{nshuffle}{epoch} = prctile(boot_percent_sig_events,[2.5 97.5],2);
                percent_multi_events1{nmethod}{nshuffle}{epoch} = boot_percent_multi_events;
                percent_multi_events_CI1{nmethod}{nshuffle}{epoch} = prctile(boot_percent_multi_events,[2.5 97.5],2);
                replay_score_compare1{nmethod}{nshuffle}{epoch} = boot_replay_score;
                replay_score_compare_CI1{nmethod}{nshuffle}{epoch} = prctile(boot_replay_score,[2.5 97.5],2);

                if nshuffle == 2
                    percent_shuffle_events1{nmethod}{epoch} = boot_percent_shuffle_events;
                    percent_shuffle_events_CI1{nmethod}{epoch} = prctile(boot_percent_shuffle_events,[2.5 97.5],2);
                end
                
                toc
            end
        end
        
    end
    log_odd_difference = log_odd_difference1;
    log_odd_difference_CI = log_odd_difference_CI1;
    percent_sig_events = percent_sig_events1;
    percent_sig_events_CI = percent_sig_events_CI1;
    percent_multi_events = percent_multi_events1;
    percent_multi_events_CI = percent_multi_events_CI1;
    percent_shuffle_events = percent_shuffle_events1;
    percent_shuffle_events_CI = percent_shuffle_events_CI1;
    replay_score_compare = replay_score_compare1;
    replay_score_compare_CI = replay_score_compare_CI1;
    
    save log_odd_difference_optimisation log_odd_difference log_odd_difference_CI percent_sig_events percent_sig_events_CI...
        percent_multi_events percent_multi_events_CI percent_shuffle_events percent_shuffle_events_CI replay_score_compare replay_score_compare_CI
    clear log_odd_difference1 log_odd_difference_CI1 percent_sig_events1 percent_sig_events_CI1 ...
        percent_multi_events1 percent_multi_events_CI1 percent_shuffle_events1 percent_shuffle_events_CI1 replay_score_compare1 replay_score_compare_CI1
    clear boot_log_odd_difference boot_percent_sig_events boot_percent_multi_events boot_percent_shuffle_events boot_replay_score
    cd ..
else
    load log_odd_difference_optimisation
    cd ..
end


cd ground_truth_original
if exist('log_odd_difference_optimisation_original.mat', 'file') ~= 2
    for nmethod = 1:length(method)
        for nshuffle = 1:length(log_pval{nmethod})
            for epoch = 1:3
                %                 multi_event_index = find(log_odd_compare{nshuffle}.multi_track == 1);
                tic
                for threshold = 1:length(p_val_threshold)
                    % original event index
                    resampled_event = epoch_index{nmethod}{nshuffle}{epoch};
                    %                     resampled_event = datasample(epoch_index{nmethod}{nshuffle}{epoch},length(epoch_index{nmethod}{nshuffle}{epoch}));

                    % Detecting events significant for track 1
                    % Find significant 'WHOLE' event
                    this_segment_event = resampled_event(find(segment_id{nmethod}{nshuffle}(1,resampled_event) == 1));
                    track_1_index = this_segment_event(find(log_pval{nmethod}{nshuffle}(1,this_segment_event) <= p_val_threshold(threshold)));
                    % Find significant half event
                    this_segment_event = resampled_event(find(segment_id{nmethod}{nshuffle}(1,resampled_event) > 1));
                    track_1_index = [track_1_index this_segment_event(find(log_pval{nmethod}{nshuffle}(1,this_segment_event) <= p_val_threshold(threshold)/2))];

                    % Detecting events significant for track 2
                    % Find significant 'WHOLE' event
                    this_segment_event = resampled_event(find(segment_id{nmethod}{nshuffle}(2,resampled_event) == 1));
                    track_2_index = this_segment_event(find(log_pval{nmethod}{nshuffle}(2,this_segment_event) <= p_val_threshold(threshold)));
                    % Find significant half event
                    this_segment_event = resampled_event(find(segment_id{nmethod}{nshuffle}(2,resampled_event) > 1));
                    track_2_index = [track_2_index this_segment_event(find(log_pval{nmethod}{nshuffle}(2,this_segment_event) <= p_val_threshold(threshold)/2))];

                    if strcmp(method{nmethod},'wcorr 1 shuffle + jump distance') % Apply jump distance threshold 0.4
                        track_1_index = track_1_index(find(max_jump{nmethod}{nshuffle}(1,track_1_index) <= 20*0.4));
                        track_2_index = track_2_index(find(max_jump{nmethod}{nshuffle}(2,track_2_index) <= 20*0.4));
                    end
                    
                    % Calculate mean replay score across all significant
                    % events
                    boot_replay_score(threshold) = mean([replay_score{nmethod}{nshuffle}(1,track_1_index) replay_score{nmethod}{nshuffle}(2,track_2_index)]);

                    % Calculate percentage of tracks significant for
                    % both tracks
                    multi_event_number = sum(ismember(track_1_index,track_2_index));
                    multi_event_percent = multi_event_number/(length(track_1_index)+length(track_2_index)-multi_event_number);

                    % Calculate mean z-scored log odds difference for all significant
                    % events (Track 1 - Track 2)
                    boot_log_odd_difference(threshold) = mean(data{nmethod}{nshuffle}(1,track_1_index)) ...
                        - mean(data{nmethod}{nshuffle}(2,track_2_index));

                    boot_percent_multi_events(threshold) = multi_event_percent;
                    % Significant event proportion (minus multitrack event number to avoid double counting)
                    boot_percent_sig_events(threshold) = (length(track_1_index)+length(track_2_index)-multi_event_number)/total_number{nshuffle}(epoch);
                    
                    if nshuffle == 2
                        % cell id shuffled mean significant event proportion
                        boot_percent_shuffle_events(threshold) = ((length(track_1_index)+length(track_2_index))/2) / total_number{nshuffle}(epoch);
                    end
                end
                
                log_odd_difference_original{nmethod}{nshuffle}{epoch} = boot_log_odd_difference;
                percent_sig_events_original{nmethod}{nshuffle}{epoch} = boot_percent_sig_events;
                percent_multi_events_original{nmethod}{nshuffle}{epoch} = boot_percent_multi_events;
                replay_score_original{nmethod}{nshuffle}{epoch} = boot_replay_score;

                if nshuffle == 2
                    percent_shuffle_events_original{nmethod}{epoch} = boot_percent_shuffle_events;
                end
                
                toc
                
            end
            
        end
    end
    
    save log_odd_difference_optimisation_original log_odd_difference_original percent_sig_events_original...
        percent_multi_events_original percent_shuffle_events_original replay_score_original
    clear boot_log_odd_difference boot_percent_sig_events boot_percent_multi_events boot_replay_score
    cd ..
    
else
    load log_odd_difference_optimisation_original
    cd ..
end

%% Multi-track proportion

p_val_threshold = round(10.^[-3:0.01:log10(0.2)],4);
p_val_threshold = flip(unique(p_val_threshold));
p_val_threshold(1) = 0.2;
low_threshold = find(ismembertol(p_val_threshold,[0.0501 0.02 0.01 0.005 0.002 0.001],eps) == 1);
Behavioural_epoches = {'PRE','RUN','POST'};
marker_shape = {'o','^','s'};
method_type = {'wcorr 1 shuffle','wcorr 1 shuffle + jump distance','wcorr 2 shuffles','wcorr 3 shuffles'...
    'spearman median spike','spearman all spikes',...
    'linear 1 shuffle','linear 2 shuffles'};
colour_line= {[127,205,187]/255,[65,182,196]/255,[34,94,168]/255,[37,52,148]/255,...
    [253,141,60]/255,[254,217,11]/255,...
    [128,0,38]/255,[227,26,28]/255};
nfig = 1;

fig = figure(2)
fig.Position = [834 116 850 885];
for epoch = 1:3
    for nmethod = 1:8
        subplot(2,2,nfig)
        y = mean(percent_multi_events{nmethod}{1}{epoch}(low_threshold,:),2);

        s(nmethod) = scatter(p_val_threshold(low_threshold),y,20,marker_shape{epoch},...
            'MarkerFaceColor',colour_line{nmethod},'MarkerEdgeColor',colour_line{nmethod},'MarkerFaceAlpha','0.5')
        hold on
        m{nmethod} = errorbar(p_val_threshold(low_threshold),y,y-percent_multi_events_CI{nmethod}{1}{epoch}(low_threshold,1),...
            percent_multi_events_CI{nmethod}{1}{epoch}(low_threshold,2)-y,'.','MarkerFaceColor',colour_line{nmethod},'MarkerEdgeColor',colour_line{nmethod})
        m{nmethod}.Color = colour_line{nmethod};
        %             hold on
        %             r(nshuffle) = rectangle('Position',pos)
        set(gca,'xscale','log')
        xticks(flip(p_val_threshold(low_threshold)))
        xticklabels({'0.001','0.002','0.005','0.01','0.02','0.05'})
        xlim([0 0.06])
        ylim([0 0.25])
        ylabel('Proportion of multi-track events')
        hold on
        xlabel('Sequenceness p value')
        title(Behavioural_epoches{epoch})
        %             title('Multitrack event proportion vs p value')
    end
    ax = gca;
    ax.FontSize = 12;
    set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
    nfig = nfig + 1;
end
legend([s(1:8)], {method_type{1:8}},'Position',[0.7 0.2 0.05 0.05])
% legend([s(5) s(6)], {method_type{5},method_type{6}},'Position',[0.7 0.2 0.05 0.05])
%     legend([s(1),s(2),s(3),s(4)], {'spike train circular shift','place field circular shift','place bin circular shift','time bin permutation'},'Position',[0.35 0.2 0.05 0.05])
%     sgtitle('Multitrack event proportion vs Sequenceness p value')
cd ground_truth_original\Figure
filename = sprintf('multi-track event proportion.pdf')
saveas(gcf,filename)
filename = sprintf('multi-track event proportion.fig')
saveas(gcf,filename)
cd ..
cd ..
clf

%% Multitrack event PRE VS RUN VS POST
p_val_threshold = round(10.^[-3:0.01:log10(0.2)],4);
p_val_threshold = flip(unique(p_val_threshold));
p_val_threshold(1) = 0.2;
low_threshold = find(ismembertol(p_val_threshold,[0.0501 0.02 0.01 0.005 0.002 0.001],eps) == 1);
Behavioural_epoches = {'PRE','RUN','POST'};
marker_shape = {'o','^','s'};
method_type = {'wcorr 1 shuffle','wcorr 1 shuffle + jump distance','wcorr 2 shuffles','wcorr 3 shuffles'...
    'spearman median spike','spearman all spikes',...
    'linear 1 shuffle','linear 2 shuffles'};

colour_line= {[0 0 0],[34,94,168]/255,[253,141,60]/255}
nfig = 1;
count = 1;
for nmethod = 1:length(method_type)
    if nfig ==1
        fig = figure
        fig.Position = [834 116 850 885];
    end

    for epoch = 1:3
        for nshuffle = 1
            subplot(2,2,nfig)
            y = mean(percent_multi_events{nmethod}{nshuffle}{epoch}(low_threshold,:),2);

            s(epoch) = scatter(p_val_threshold(low_threshold),y,20,marker_shape{epoch},...
                'MarkerFaceColor',colour_line{epoch},'MarkerEdgeColor',colour_line{epoch},'MarkerFaceAlpha','0.5')
            hold on
            m{epoch} = errorbar(p_val_threshold(low_threshold),y,y-percent_multi_events_CI{nmethod}{1}{epoch}(low_threshold,1),...
                percent_multi_events_CI{nmethod}{1}{epoch}(low_threshold,2)-y,'.','MarkerFaceColor',colour_line{epoch},'MarkerEdgeColor',colour_line{epoch})
            m{epoch}.Color = colour_line{epoch};
            %             hold on
            %             r(nshuffle) = rectangle('Position',pos)

            set(gca,'xscale','log')
            xticks(flip(p_val_threshold(low_threshold)))
            xticklabels({'0.001','0.002','0.005','0.01','0.02','0.05'})
            xlim([0 0.06])
            ylim([0 0.25])
            ylabel('Proportion of multi-track events')
            hold on
            xlabel('Sequenceness p value')

            title(method_type{nmethod})
        end
    end
    ax = gca;
    ax.FontSize = 12;
    set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
    nfig = nfig + 1;
    legend([s(1:3)], {Behavioural_epoches{1:3}},'Position',[0.7 0.2 0.05 0.05])


    for epoch = 1:3
        for nshuffle = 1
            subplot(2,2,nfig)
            y = mean(percent_shuffle_events{nmethod}{epoch}(low_threshold,:),2);

            s(epoch) = scatter(p_val_threshold(low_threshold),y,20,marker_shape{epoch},...
                'MarkerFaceColor',colour_line{epoch},'MarkerEdgeColor',colour_line{epoch},'MarkerFaceAlpha','0.5')
            hold on
            m{epoch} = errorbar(p_val_threshold(low_threshold),y,y-percent_shuffle_events_CI{nmethod}{epoch}(low_threshold,1),...
                percent_shuffle_events_CI{nmethod}{epoch}(low_threshold,2)-y,'.','MarkerFaceColor',colour_line{epoch},'MarkerEdgeColor',colour_line{epoch})
            m{epoch}.Color = colour_line{epoch};
            %             hold on
            %             r(nshuffle) = rectangle('Position',pos)

            set(gca,'xscale','log')
            xticks(flip(p_val_threshold(low_threshold)))
            xticklabels({'0.001','0.002','0.005','0.01','0.02','0.05'})
            xlim([0 0.06])
            ylim([0 0.25])
            ylabel('Mean false positive rates')
            hold on
            xlabel('Sequenceness p value')

            title(method_type{nmethod})
        end
    end
    ax = gca;
    ax.FontSize = 12;
    set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
    nfig = nfig + 1;
    legend([s(1:3)], {Behavioural_epoches{1:3}},'Position',[0.7 0.2 0.05 0.05])




    if nfig==5
        nfig = 1;
        cd ground_truth_original\Figure
        filename = sprintf('multi-track event proportion PRE vs POST %i.pdf',count)
        saveas(gcf,filename)
        filename = sprintf('multi-track event proportio PRE vs POST %i.fig',count)
        saveas(gcf,filename)
        cd ..
        cd ..
        count = count + 1;
%         clf
    end
end
%% Multitrack event vs false positive rate
p_val_threshold = round(10.^[-3:0.01:log10(0.2)],4);
p_val_threshold = flip(unique(p_val_threshold));
p_val_threshold(1) = 0.2;
low_threshold = find(ismembertol(p_val_threshold,[0.0501 0.02 0.01 0.005 0.002 0.001],eps) == 1);
Behavioural_epoches = {'PRE','RUN','POST'};
marker_shape = {'o','^','s'};
method_type = {'wcorr 1 shuffle','wcorr 1 shuffle + jump distance','wcorr 2 shuffles','wcorr 3 shuffles'...
    'spearman median spike','spearman all spikes',...
    'linear 1 shuffle','linear 2 shuffles'};

colour_line= {[127,205,187]/255,[65,182,196]/255,[34,94,168]/255,[37,52,148]/255,...
    [253,141,60]/255,[254,217,11]/255,...
    [128,0,38]/255,[227,26,28]/255};
nfig = 1;
count = 1;
clear text
fig = figure
fig.Position = [834 116 850 885];
subplot(2,2,1)
for epoch = 1:3
    for nmethod = 1:length(method)
        x_CI = percent_shuffle_events_CI{nmethod}{epoch}(low_threshold(1),:);
        x = mean(percent_shuffle_events{nmethod}{epoch}(low_threshold(1),:),2);

        y_CI = percent_multi_events_CI{nmethod}{1}{epoch}(low_threshold(1),:);
        y = mean(percent_multi_events{nmethod}{1}{epoch}(low_threshold(1),:),2);

        s(count) = scatter(x,y,20,marker_shape{epoch},...
            'MarkerFaceColor',colour_line{nmethod},'MarkerEdgeColor',colour_line{nmethod},'MarkerFaceAlpha','0.5')
        %             hold on
        %             m{nshuffle} = errorbar(p_val_threshold(low_threshold),y,y-percent_multi_events_CI{nmethod}{nshuffle}{epoch}(low_threshold,1),...
        %                 percent_multi_events_CI{nmethod}{nshuffle}{epoch}(low_threshold,2)-y,'.','MarkerFaceColor',colour_line{nshuffle},'MarkerEdgeColor',colour_line{nshuffle})
        %             m{nshuffle}.Color = colour_line{nshuffle};
        hold on
        error_x = [x_CI(1) x_CI(2) x_CI(2) x_CI(1)];
        error_y = [y_CI(1) y_CI(1) y_CI(2) y_CI(2)];

        patch(error_x,error_y,colour_line{nmethod},'EdgeColor',colour_line{nmethod},'FaceAlpha','0.2','EdgeAlpha','0.2')

        xlabel('Mean false positive rates')
        hold on
        ylabel('Mean proportion of multi-track events')
        %         title(Behavioural_epoches{epoch})
        %         title('Original p value =< 0.05')

        text{count} = sprintf('%s %s'...
            ,method_type{nmethod},Behavioural_epoches{epoch});
        count = count + 1;

    end
    ax = gca;
    ax.FontSize = 12;
    set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
    ylim([0 0.25])
    xlim([0 0.25])
end
plot([0 0.25],[0 0.25],'k--')
sgtitle('Alpha level = 0.05')
% legend([s(:)], {method{:}},'Position',[0.7 0.2 0.05 0.05])
legend([s(1:24)],...
    {text{1:24}},'Position',[0.7 0.3 0.05 0.05])

cd ground_truth_original\Figure
filename = 'multi-track event proportion vs FPR at alpha 0.05.pdf';
saveas(gcf,filename)
filename = 'multi-track event proportion vs FPR at alpha 0.05.fig'
saveas(gcf,filename)
cd ..
cd ..

% Mean FPR 0.05
nfig = 1;
count = 1;
clear text
fig = figure
fig.Position = [834 116 850 885];
for epoch = 1:3
    for nmethod = 1:length(method)
        [c index(nmethod)] = min(abs(mean(percent_shuffle_events{nmethod}{epoch},2) - 0.05));
        x_CI = percent_shuffle_events_CI{nmethod}{epoch}(index(nmethod),:);
        x = mean(percent_shuffle_events{nmethod}{epoch}(index(nmethod),:),2);

        y_CI = percent_multi_events_CI{nmethod}{1}{epoch}(index(nmethod),:);
        y = mean(percent_multi_events{nmethod}{1}{epoch}(index(nmethod),:),2);

        s(count) = scatter(x,y,20,marker_shape{epoch},...
            'MarkerFaceColor',colour_line{nmethod},'MarkerEdgeColor',colour_line{nmethod},'MarkerFaceAlpha','0.5')
        %             hold on
        %             m{nshuffle} = errorbar(p_val_threshold(low_threshold),y,y-percent_multi_events_CI{nmethod}{nshuffle}{epoch}(low_threshold,1),...
        %                 percent_multi_events_CI{nmethod}{nshuffle}{epoch}(low_threshold,2)-y,'.','MarkerFaceColor',colour_line{nshuffle},'MarkerEdgeColor',colour_line{nshuffle})
        %             m{nshuffle}.Color = colour_line{nshuffle};
        hold on
        error_x = [x_CI(1) x_CI(2) x_CI(2) x_CI(1)];
        error_y = [y_CI(1) y_CI(1) y_CI(2) y_CI(2)];

%         patch(error_x,error_y,colour_line{nmethod},'EdgeColor',colour_line{nmethod},'FaceAlpha','0.2','EdgeAlpha','0.2')

        xlabel('Mean false positive rates')
        hold on
        ylabel('Mean proportion of multi-track events')
        %         title(Behavioural_epoches{epoch})
        %         title('Original p value =< 0.05')

        text{count} = sprintf('%s %s'...
            ,method_type{nmethod},Behavioural_epoches{epoch});
        count = count + 1;

    end
    ax = gca;
    ax.FontSize = 12;
    set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
    ylim([0 0.25])
    xlim([0 0.25])
end
sgtitle('Mean false positive rate = 0.05')
% legend([s(:)], {method{:}},'Position',[0.7 0.2 0.05 0.05])
legend([s(1:24)],...
    {text{1:24}})

cd ground_truth_original\Figure
filename = 'multi-track event proportion vs FPR at matching FPR.pdf';
saveas(gcf,filename)
filename = 'multi-track event proportion vs FPR at matching FPR.fig'
saveas(gcf,filename)
cd ..
cd ..

%% Multitrack event original vs cell id randomised dataset
p_val_threshold = round(10.^[-3:0.01:log10(0.2)],4);
p_val_threshold = flip(unique(p_val_threshold));
p_val_threshold(1) = 0.2;
low_threshold = find(ismembertol(p_val_threshold,[0.0501 0.02 0.01 0.005 0.002 0.001],eps) == 1);
Behavioural_epoches = {'PRE','RUN','POST'};
marker_shape = {'o','^','s'};
method_type = {'wcorr 1 shuffle','wcorr 1 shuffle + jump distance','wcorr 2 shuffles','wcorr 3 shuffles'...
    'spearman median spike','spearman all spikes',...
    'linear 1 shuffle','linear 2 shuffles'};

colour_line= {[127,205,187]/255,[65,182,196]/255,[34,94,168]/255,[37,52,148]/255,...
    [253,141,60]/255,[254,217,11]/255,...
    [128,0,38]/255,[227,26,28]/255};
nfig = 1;

for nmethod = 1:length(method_type)
    fig = figure(nmethod)
fig.Position = [834 116 850 885];
    for epoch = 1:3
        for nshuffle = 1:2
            subplot(2,2,epoch)
            y = mean(percent_multi_events{nmethod}{nshuffle}{epoch}(low_threshold,:),2);
            if nshuffle == 1
                s(nshuffle) = scatter(p_val_threshold(low_threshold),y,20,marker_shape{epoch},...
                    'MarkerFaceColor',colour_line{nmethod},'MarkerEdgeColor',colour_line{nmethod},'MarkerFaceAlpha','0.5')
                hold on
                m{nshuffle} = errorbar(p_val_threshold(low_threshold),y,y-percent_multi_events_CI{nmethod}{1}{epoch}(low_threshold,1),...
                    percent_multi_events_CI{nmethod}{1}{epoch}(low_threshold,2)-y,'.','MarkerFaceColor',colour_line{nmethod},'MarkerEdgeColor',colour_line{nmethod})
                m{nshuffle}.Color = colour_line{nmethod};
                %             hold on
                %             r(nshuffle) = rectangle('Position',pos)
            else

                s(nshuffle) = scatter(p_val_threshold(low_threshold),y,20,marker_shape{epoch},...
                    'MarkerFaceColor',colour_line{nshuffle},'MarkerEdgeColor','k','MarkerFaceAlpha','0.5')
                hold on
                m{nshuffle} = errorbar(p_val_threshold(low_threshold),y,y-percent_multi_events_CI{nmethod}{2}{epoch}(low_threshold,1),...
                    percent_multi_events_CI{nmethod}{2}{epoch}(low_threshold,2)-y,'.','MarkerFaceColor','k','MarkerEdgeColor','k')
                m{nshuffle}.Color = 'k';
                %             hold on
                %             r(nshuffle) = rectangle('Position',pos)
            end
            set(gca,'xscale','log')
            xticks(flip(p_val_threshold(low_threshold)))
            xticklabels({'0.001','0.002','0.005','0.01','0.02','0.05'})
            xlim([0 0.06])
            ylim([0 0.25])
            ylabel('Proportion of multi-track events')
            hold on
            xlabel('Sequenceness p value')
            title(Behavioural_epoches{epoch})
            %             title('Multitrack event proportion vs p value')
        end
    end
    ax = gca;
    ax.FontSize = 12;
    set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
    nfig = nfig + 1;
    legend([s(1:2)], {'Original','Cell ID randomized dataset'},'Position',[0.7 0.2 0.05 0.05])
    sgtitle(method_type{nmethod})
end