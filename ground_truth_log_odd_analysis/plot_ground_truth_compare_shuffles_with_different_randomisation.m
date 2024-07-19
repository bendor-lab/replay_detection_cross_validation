function [] = plot_ground_truth_compare_shuffles_with_different_randomisation(folders,option,method)

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

total_number{3} = total_number{2};
total_number{4} = total_number{2};
total_number{5} = total_number{2};

index = [];
epoch_index = [];
for nmethod = 1:length(method)
    cd ground_truth_original
    if strcmp(method{nmethod},'wcorr 1 shuffle')
        load log_odd_wcorr_place_bin_circular_shift
        log_odd_compare{nmethod}{1} = log_odd;
        load log_odd_wcorr_place_bin_circular_shift_global_remapped
        log_odd_compare{nmethod}{2} = log_odd;
        load log_odd_wcorr_place_bin_circular_shift_cross_experiment_shuffled
        log_odd_compare{nmethod}{3} = log_odd;
        load log_odd_wcorr_place_bin_circular_shift_place_field_shifted
        log_odd_compare{nmethod}{4} = log_odd;
        load log_odd_wcorr_place_bin_circular_shift_spike_train_shifted
        log_odd_compare{nmethod}{5} = log_odd;


    elseif strcmp(method{nmethod},'wcorr 1 shuffle + jump distance')
        load log_odd_wcorr_place_bin_circular_shift
        log_odd_compare{nmethod}{1} = log_odd;
        load log_odd_wcorr_place_bin_circular_shift_global_remapped
        log_odd_compare{nmethod}{2} = log_odd;
        load log_odd_wcorr_place_bin_circular_shift_cross_experiment_shuffled
        log_odd_compare{nmethod}{3} = log_odd;
        load log_odd_wcorr_place_bin_circular_shift_place_field_shifted
        log_odd_compare{nmethod}{4} = log_odd;
        load log_odd_wcorr_place_bin_circular_shift_spike_train_shifted
        log_odd_compare{nmethod}{5} = log_odd;

    elseif strcmp(method{nmethod},'wcorr 2 shuffles')
        load log_odd_wcorr_PRE_place_POST_time
        log_odd_compare{nmethod}{1} = log_odd;
        load log_odd_wcorr_PRE_place_POST_time_global_remapped
        log_odd_compare{nmethod}{2} = log_odd;
        load log_odd_wcorr_PRE_place_POST_time_cross_experiment_shuffled
        log_odd_compare{nmethod}{3} = log_odd;
        load log_odd_wcorr_PRE_place_POST_time_place_field_shifted
        log_odd_compare{nmethod}{4} = log_odd;
        load log_odd_wcorr_PRE_place_POST_time_spike_train_shifted
        log_odd_compare{nmethod}{5} = log_odd;

    elseif strcmp(method{nmethod},'wcorr 3 shuffles')
        load log_odd_wcorr
        log_odd_compare{nmethod}{1} = log_odd;
        load log_odd_wcorr_global_remapped
        log_odd_compare{nmethod}{2} = log_odd;
        load log_odd_wcorr_cross_experiment_shuffled
        log_odd_compare{nmethod}{3} = log_odd;
        load log_odd_wcorr_place_field_shifted
        log_odd_compare{nmethod}{4} = log_odd;
        load log_odd_wcorr_spike_train_shifted
        log_odd_compare{nmethod}{5} = log_odd;

    elseif strcmp(method{nmethod},'spearman median spike')
        load log_odd_spearman
        log_odd_compare{nmethod}{1} = log_odd;
        load log_odd_spearman_global_remapped
        log_odd_compare{nmethod}{2} = log_odd;
        load log_odd_spearman_cross_experiment_shuffled
        log_odd_compare{nmethod}{3} = log_odd;
        load log_odd_spearman_place_field_shifted
        log_odd_compare{nmethod}{4} = log_odd;
        load log_odd_spearman_spike_train_shifted
        log_odd_compare{nmethod}{5} = log_odd;

    elseif strcmp(method{nmethod},'spearman all spikes')
        load log_odd_spearman_all_spikes
        log_odd_compare{nmethod}{1} = log_odd;
        load log_odd_spearman_all_spikes_global_remapped
        log_odd_compare{nmethod}{2} = log_odd;
        load log_odd_spearman_all_spikes_cross_experiment_shuffled
        log_odd_compare{nmethod}{3} = log_odd;
        load log_odd_spearman_all_spikes_place_field_shifted
        log_odd_compare{nmethod}{4} = log_odd;
        load log_odd_spearman_all_spikes_spike_train_shifted
        log_odd_compare{nmethod}{5} = log_odd;

    elseif strcmp(method{nmethod},'linear 1 shuffle')
        load log_odd_linear_place_bin_circular_shift
        log_odd_compare{nmethod}{1} = log_odd;
        load log_odd_linear_place_bin_circular_shift_global_remapped
        log_odd_compare{nmethod}{2} = log_odd;
        load log_odd_linear_place_bin_circular_shift_cross_experiment_shuffled
        log_odd_compare{nmethod}{3} = log_odd;
        load log_odd_linear_place_bin_circular_shift_place_field_shifted
        log_odd_compare{nmethod}{4} = log_odd;
        load log_odd_linear_place_bin_circular_shift_spike_train_shifted
        log_odd_compare{nmethod}{5} = log_odd;

    elseif strcmp(method{nmethod},'linear 2 shuffles')
        load log_odd_linear_PRE_place_POST_time
        log_odd_compare{nmethod}{1} = log_odd;
        load log_odd_linear_PRE_place_POST_time_global_remapped
        log_odd_compare{nmethod}{2} = log_odd;
        load log_odd_linear_PRE_place_POST_time_cross_experiment_shuffled
        log_odd_compare{nmethod}{3} = log_odd;
        load log_odd_linear_PRE_place_POST_time_place_field_shifted
        log_odd_compare{nmethod}{4} = log_odd;
        load log_odd_linear_PRE_place_POST_time_spike_train_shifted
        log_odd_compare{nmethod}{5} = log_odd;

    end
    cd ..

    % Get index for PRE,  RUN, POST
    states = [-1 0 1 2]; % PRE, POST, RUN Track 1, RUN Track 2
    
    for nshuffle = 1:length(log_odd_compare{nmethod})
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
        data{nmethod}{3} = log_odd_compare{nmethod}{3}.common_zscored.cross_experiment_shuffled_original;
        data{nmethod}{4} = log_odd_compare{nmethod}{4}.common_zscored.place_field_shifted_original;
        data{nmethod}{5} = log_odd_compare{nmethod}{5}.common_zscored.spike_train_shifted_original;

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
if exist('log_odd_difference_compare_different_randomisation.mat', 'file') ~= 2
  
    for nmethod = 1:length(method)
        for epoch = 1:3
            for nshuffle = 1:length(log_pval{nmethod})
                tic
                
                for threshold = 1:length(p_val_threshold)
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

                        if  nshuffle ~= 1
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

                if nshuffle ~= 1
                    percent_shuffle_events1{nmethod}{epoch}{nshuffle} = boot_percent_shuffle_events;
                    percent_shuffle_events_CI1{nmethod}{epoch}{nshuffle} = prctile(boot_percent_shuffle_events,[2.5 97.5],2);
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
    
    save log_odd_difference_compare_different_randomisation log_odd_difference log_odd_difference_CI percent_sig_events percent_sig_events_CI...
        percent_multi_events percent_multi_events_CI percent_shuffle_events percent_shuffle_events_CI replay_score_compare replay_score_compare_CI
    clear log_odd_difference1 log_odd_difference_CI1 percent_sig_events1 percent_sig_events_CI1 ...
        percent_multi_events1 percent_multi_events_CI1 percent_shuffle_events1 percent_shuffle_events_CI1 replay_score_compare1 replay_score_compare_CI1
    clear boot_log_odd_difference boot_percent_sig_events boot_percent_multi_events boot_percent_shuffle_events boot_replay_score
    cd ..
else
    load log_odd_difference_compare_different_randomisation
    cd ..
end


cd ground_truth_original
if exist('log_odd_difference_compare_different_randomisation_original.mat', 'file') ~= 2
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
                    
                    if nshuffle ~= 2
                        % cell id shuffled mean significant event proportion
                        boot_percent_shuffle_events(threshold) = ((length(track_1_index)+length(track_2_index))/2) / total_number{nshuffle}(epoch);
                    end
                end
                
                log_odd_difference_original{nmethod}{nshuffle}{epoch} = boot_log_odd_difference;
                percent_sig_events_original{nmethod}{nshuffle}{epoch} = boot_percent_sig_events;
                percent_multi_events_original{nmethod}{nshuffle}{epoch} = boot_percent_multi_events;
                replay_score_original{nmethod}{nshuffle}{epoch} = boot_replay_score;

                if nshuffle ~= 1
                    percent_shuffle_events_original{nmethod}{epoch}{nshuffle} = boot_percent_shuffle_events;
                end
                
                toc
                
            end
            
        end
    end
    
    save log_odd_difference_compare_different_randomisation_original log_odd_difference_original percent_sig_events_original...
        percent_multi_events_original percent_shuffle_events_original replay_score_original
    clear boot_log_odd_difference boot_percent_sig_events boot_percent_multi_events boot_replay_score
    cd ..
    
else
    load log_odd_difference_compare_different_randomisation_original
    cd ..
end


%% Spearman optimisation


% Inset for false positive rate at p value 0.05 and p value for false
% positive rate 0.05

fig = figure(3)
fig.Position = [834 116 575 531];
count = 1;

subplot(2,2,1)
for nmethod = [5 6]
    for epoch = 2:3
        mean_false_positive_rate(count) = mean(percent_shuffle_events{nmethod}{epoch}(low_threshold(1),:),2);
        %         b = prctile(percent_shuffle_events_POST_ripple{nmethod}{epoch}(low_threshold(1),:),[2.5 97.5])

        s(nmethod) = scatter(epoch-1.5,mean_false_positive_rate(count),20,marker_shape{epoch},...
            'MarkerFaceColor',colour_line{nmethod},'MarkerEdgeColor',colour_line{nmethod},'MarkerFaceAlpha','0.5')
        hold on
        m{nmethod} = errorbar(epoch-1.5,mean_false_positive_rate(count),mean_false_positive_rate(count)-percent_shuffle_events_CI{nmethod}{epoch}(low_threshold(1),1),...
            percent_shuffle_events_CI{nmethod}{epoch}(low_threshold(1),2)-mean_false_positive_rate(count),'.','MarkerFaceColor',colour_line{nmethod},'MarkerEdgeColor',colour_line{nmethod})
        m{nmethod}.Color = colour_line{nmethod};
        %         set(gca,'yscale','log')
        count = count + 1;
    end
    ax = gca;
    ax.FontSize = 12;
    set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
end
ylabel('mean false positive rate')
pbaspect([1 3 1])
ylim([0 0.2])
xlim([0 2])
set(gca,'xtick',[]);

subplot(2,2,2)
count = 1;
for nmethod = [5 6]
    for epoch = 2:3
        [xx index] = min(abs(mean(percent_shuffle_events{nmethod}{epoch}{ndataset},2) - 0.05));
        adjusted_pvalue(count) = p_val_threshold(index);

        s(nmethod) = scatter(epoch-1.5,adjusted_pvalue(count),20,marker_shape{epoch},...
            'MarkerFaceColor',colour_line{nmethod},'MarkerEdgeColor',colour_line{nmethod},'MarkerFaceAlpha','0.5')
        hold on
        %         m{nmethod} = errorbar(p_val_threshold(low_threshold),y,y-percent_shuffle_events_POST_ripple_CI{nmethod}{epoch}(low_threshold,1),...
        %             percent_shuffle_events_POST_ripple_CI{nmethod}{epoch}(low_threshold,2)-y,'.','MarkerFaceColor',colour_line{nmethod},'MarkerEdgeColor',colour_line{nmethod})
        %         m{nmethod}.Color = colour_line{nmethod};
        %         set(gca,'yscale','log')
        count = count + 1;
    end
    ax = gca;
    ax.FontSize = 12;
    set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
end
ylabel('Adjusted p value')
pbaspect([1 3 1])
ylim([0 0.05])
xlim([0 2])
set(gca,'xtick',[]);
cd ground_truth_original\Figure
filename = sprintf('spearman optimisation false positive rate and p value inset (global remapped).pdf')
saveas(gcf,filename)
filename = sprintf('spearman optimisation false positive rate and p value inset (POST ripple).fig')
saveas(gcf,filename)
cd ..
cd ..
clf




% Sig proportion and log odd at p value 0.05
fig = figure(3)
fig.Position = [834 116 850 885];
subplot(2,2,1)
count = 1;
for nmethod = [5 6]
    for epoch = 2:3
        x = mean(mean(log_odd_difference{nmethod}{1}{epoch}(low_threshold(1),:),2),2);
        y = mean(percent_sig_events{nmethod}{1}{epoch}(low_threshold(1),:),2);
        s(count) = scatter(x,y,20,marker_shape{epoch},...
            'MarkerFaceColor',colour_line{nmethod},'MarkerEdgeColor',colour_line{nmethod},'MarkerFaceAlpha','0.5')
        %             hold on
        %             m{nshuffle} = errorbar(p_val_threshold(low_threshold),y,y-percent_multi_events_CI{nmethod}{nshuffle}{epoch}(low_threshold,1),...
        %                 percent_multi_events_CI{nmethod}{nshuffle}{epoch}(low_threshold,2)-y,'.','MarkerFaceColor',colour_line{nshuffle},'MarkerEdgeColor',colour_line{nshuffle})
        %             m{nshuffle}.Color = colour_line{nshuffle};
        hold on
        error_x = [log_odd_difference_CI{nmethod}{1}{epoch}(low_threshold(1),1) log_odd_difference_CI{nmethod}{1}{epoch}(low_threshold(1),2)...
            log_odd_difference_CI{nmethod}{1}{epoch}(low_threshold(1),2) log_odd_difference_CI{nmethod}{1}{epoch}(low_threshold(1),1)]
        error_y = [percent_sig_events_CI{nmethod}{1}{epoch}(low_threshold(1),1) percent_sig_events_CI{nmethod}{1}{epoch}(low_threshold(1),1)...
            percent_sig_events_CI{nmethod}{1}{epoch}(low_threshold(1),2) percent_sig_events_CI{nmethod}{1}{epoch}(low_threshold(1),2)]

        patch(error_x,error_y,colour_line{nmethod},'EdgeColor',colour_line{nmethod},'FaceAlpha','0.2','EdgeAlpha','0.2')
        text{count} = sprintf('%s p =< %.3f (proportion = %.3f)',method_type{nmethod},p_val_threshold(low_threshold(1)),mean(percent_shuffle_events{nmethod}{epoch}(low_threshold(1),:)));

        ylabel('Proportion of significant events')
        hold on
        xlabel('mean log odds difference')
        %             title(Behavioural_epoches{epoch})
        title('Original p value =< 0.05')

        count = count + 1;
    end
        ax = gca;
    ax.FontSize = 12;
    set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
end
legend([s(1:4)], {text{1:4}})

% xlim([-1 3])
ylim([0 0.56])

% Sig proportion and log odd at shuffle-corrected p value 0.05
subplot(2,2,2)
count = 1;
%         tempt = find(mean(percent_multi_events{nmethod}{nshuffle}{epoch},2) <= 0.05);
for nmethod = [5 6]
    for epoch = 2:3
        [c index(nmethod)] = min(abs(mean(percent_shuffle_events{nmethod}{epoch},2) - 0.05));
%          [c index(nmethod)] = min(abs(mean(percent_multi_events{nmethod}{1}{epoch},2) - 0.05));
        x = mean(mean(log_odd_difference{nmethod}{1}{epoch}(index(nmethod),:),2),2);
        y = mean(percent_sig_events{nmethod}{1}{epoch}(index(nmethod),:),2);
        s(count) = scatter(x,y,20,marker_shape{epoch},...
            'MarkerFaceColor',colour_line{nmethod},'MarkerEdgeColor',colour_line{nmethod},'MarkerFaceAlpha','0.5')
        %             hold on
        %             m{nshuffle} = errorbar(p_val_threshold(low_threshold),y,y-percent_multi_events_CI{nmethod}{nshuffle}{epoch}(low_threshold,1),...
        %                 percent_multi_events_CI{nmethod}{nshuffle}{epoch}(low_threshold,2)-y,'.','MarkerFaceColor',colour_line{nshuffle},'MarkerEdgeColor',colour_line{nshuffle})
        %             m{nshuffle}.Color = colour_line{nshuffle};
        hold on
        error_x = [log_odd_difference_CI{nmethod}{1}{epoch}(index(nmethod),1) log_odd_difference_CI{nmethod}{1}{epoch}(index(nmethod),2)...
            log_odd_difference_CI{nmethod}{1}{epoch}(index(nmethod),2) log_odd_difference_CI{nmethod}{1}{epoch}(index(nmethod),1)]
        error_y = [percent_sig_events_CI{nmethod}{1}{epoch}(index(nmethod),1) percent_sig_events_CI{nmethod}{1}{epoch}(index(nmethod),1)...
            percent_sig_events_CI{nmethod}{1}{epoch}(index(nmethod),2) percent_sig_events_CI{nmethod}{1}{epoch}(index(nmethod),2)]

        patch(error_x,error_y,colour_line{nmethod},'EdgeColor',colour_line{nmethod},'FaceAlpha','0.2','EdgeAlpha','0.2')
        text{count} = sprintf('%s p =< %.3f (proportion = %.3f)',method_type{nmethod},p_val_threshold(index(nmethod)),mean(percent_shuffle_events{nmethod}{epoch}(index(nmethod),:)));
%         text{nmethod} = sprintf('%s p =< %.3f (proportion = %.3f)',shuffle_type{nmethod},p_val_threshold(index(nmethod)),mean(percent_multi_events{nmethod}{1}{epoch}(index(nmethod),:)));

        ylabel('Proportion of significant events')
        hold on
        xlabel('mean log odd difference')
        %             title(Behavioural_epoches{epoch})
        title('Equivalent p value when cell-id shuffled sig proportion = 0.05')
        count = count+ 1;
    end
    ax = gca;
    ax.FontSize = 12;
    set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
end
%         plot([0 0.1],[0 0.3],'r')
legend([s(1:4)], {text{1:4}})
% xlim([-1 3])
ylim([0 0.56])


cd ground_truth_original\Figure
filename = sprintf('spearman optimisation original vs shuffle-corrected p value.pdf')
saveas(gcf,filename)
filename = sprintf('spearman optimisation original vs shuffle-corrected p value.fig')
saveas(gcf,filename)
cd ..
cd ..
clf


%% Linear optimisation

p_val_threshold = round(10.^[-3:0.01:log10(0.2)],4);
% p_val_threshold = 0.001:0.0001:0.2;
p_val_threshold = flip(unique(p_val_threshold));
p_val_threshold(1) = 0.2;
low_threshold = find(ismembertol(p_val_threshold,[0.0501 0.02 0.01 0.005 0.002 0.001],eps) == 1);

colour_line= {[127,205,187]/255,[65,182,196]/255,[34,94,168]/255,[37,52,148]/255,...
    [253,141,60]/255,[254,217,11]/255,...
    [128,0,38]/255,[227,26,28]/255};

alpha_level = linspace(0.2,0.7,length(low_threshold));
fig = figure(2)
fig.Position = [834 116 850 700];
nfig = 1;

for epoch = 2:3
    subplot(2,2,nfig)
    for nmethod = [7 8]
        x = mean(log_odd_difference{nmethod}{1}{epoch},2)';
        y = mean(percent_sig_events{nmethod}{1}{epoch},2)';

        for threshold = 1:length(low_threshold)
            scatter(x(low_threshold(threshold)),y(low_threshold(threshold)),'filled','MarkerFaceColor',colour_line{nmethod},'MarkerFaceAlpha',alpha_level(threshold),'MarkerEdgeColor',colour_line{nmethod})
            hold on
        end

        UCI = log_odd_difference_CI{nmethod}{1}{epoch}(:,2)';
        UCI(isnan(x)) = [];
        LCI = log_odd_difference_CI{nmethod}{1}{epoch}(:,1)';
        LCI(isnan(x)) = [];
        y(isnan(x)) = [];
        x(isnan(x)) = [];
        p(nmethod) = patch([LCI fliplr(UCI)], [y fliplr(y)], colour_line{nmethod},'FaceAlpha','0.2','LineStyle','none');
        hold on

        xlabel('mean log odd difference')
        ylabel('Proportion of all replay events')
    end

    ylim([0 1])
    title(Behavioural_epoches{epoch});
    nfig = nfig + 1;
end
%     legend([p(1) p(3) p(6)], {'Shuffle 3','Shuffle 2+3','Shuffle 1+2+3'})
legend([p(7),p(8)], {method_type{7},method_type{8}},'Position',[0.7 0.3 0.05 0.05])

cd ground_truth_original\Figure
filename = sprintf('linear optimisation comparisions CI.pdf')
saveas(gcf,filename)
filename = sprintf('linear optimisation comparisions CI.fig')
saveas(gcf,filename)
cd ..
cd ..
clf
%



p_val_threshold = round(10.^[-3:0.01:log10(0.2)],4);
p_val_threshold = flip(unique(p_val_threshold));
p_val_threshold(1) = 0.2;
low_threshold = find(ismembertol(p_val_threshold,[0.0501 0.02 0.01 0.005 0.002 0.001],eps) == 1);
Behavioural_epoches = {'PRE','RUN','POST'};
marker_shape = {'o','^','s'};
nfig = 1;

fig = figure(2)
fig.Position = [834 116 850 885];
for epoch = 2:3
    for nmethod = [7 8]
        subplot(2,2,nfig)
        y = mean(percent_shuffle_events{nmethod}{epoch}(low_threshold,:),2);

        s(nmethod) = scatter(p_val_threshold(low_threshold),y,20,marker_shape{epoch},...
            'MarkerFaceColor',colour_line{nmethod},'MarkerEdgeColor',colour_line{nmethod},'MarkerFaceAlpha','0.5')
        hold on
        m{nmethod} = errorbar(p_val_threshold(low_threshold),y,y-percent_shuffle_events_CI{nmethod}{epoch}(low_threshold,1),...
            percent_shuffle_events_CI{nmethod}{epoch}(low_threshold,2)-y,'.','MarkerFaceColor',colour_line{nmethod},'MarkerEdgeColor',colour_line{nmethod})
        m{nmethod}.Color = colour_line{nmethod};
        %             hold on
        %             r(nshuffle) = rectangle('Position',pos)
        set(gca,'xscale','log')
        xticks(flip(p_val_threshold(low_threshold)))
        xticklabels({'0.001','0.002','0.005','0.01','0.02','0.05'})
        xlim([0 0.06])
        ylim([0 0.2])
        ylabel('Proportion of cell-id shuffled events')
        hold on
        xlabel('Sequenceness p value')
        title(Behavioural_epoches{epoch})
        %             title('Multitrack event proportion vs p value')
    end
    nfig = nfig + 1;
end
legend([s(7) s(8)], {method_type{7},method_type{8}},'Position',[0.7 0.2 0.05 0.05])
%     legend([s(1),s(2),s(3),s(4)], {'spike train circular shift','place field circular shift','place bin circular shift','time bin permutation'},'Position',[0.35 0.2 0.05 0.05])
%     sgtitle('Multitrack event proportion vs Sequenceness p value')
cd ground_truth_original\Figure
filename = sprintf('linear optimisation cell-id shuffled sig proportion.pdf')
saveas(gcf,filename)
filename = sprintf('linear optimisation cell-id shuffled sig proportion.fig')
saveas(gcf,filename)
cd ..
cd ..
clf



% Sig proportion and log odd at p value 0.05
fig = figure(3)
fig.Position = [834 116 850 885];
subplot(2,2,1)
count = 1;
for nmethod = [7 8]
    for epoch = 2:3
        x = mean(mean(log_odd_difference{nmethod}{1}{epoch}(low_threshold(1),:),2),2);
        y = mean(percent_sig_events{nmethod}{1}{epoch}(low_threshold(1),:),2);
        s(count) = scatter(x,y,20,marker_shape{epoch},...
            'MarkerFaceColor',colour_line{nmethod},'MarkerEdgeColor',colour_line{nmethod},'MarkerFaceAlpha','0.5')
        %             hold on
        %             m{nshuffle} = errorbar(p_val_threshold(low_threshold),y,y-percent_multi_events_CI{nmethod}{nshuffle}{epoch}(low_threshold,1),...
        %                 percent_multi_events_CI{nmethod}{nshuffle}{epoch}(low_threshold,2)-y,'.','MarkerFaceColor',colour_line{nshuffle},'MarkerEdgeColor',colour_line{nshuffle})
        %             m{nshuffle}.Color = colour_line{nshuffle};
        hold on
        error_x = [log_odd_difference_CI{nmethod}{1}{epoch}(low_threshold(1),1) log_odd_difference_CI{nmethod}{1}{epoch}(low_threshold(1),2)...
            log_odd_difference_CI{nmethod}{1}{epoch}(low_threshold(1),2) log_odd_difference_CI{nmethod}{1}{epoch}(low_threshold(1),1)]
        error_y = [percent_sig_events_CI{nmethod}{1}{epoch}(low_threshold(1),1) percent_sig_events_CI{nmethod}{1}{epoch}(low_threshold(1),1)...
            percent_sig_events_CI{nmethod}{1}{epoch}(low_threshold(1),2) percent_sig_events_CI{nmethod}{1}{epoch}(low_threshold(1),2)]

        patch(error_x,error_y,colour_line{nmethod},'EdgeColor',colour_line{nmethod},'FaceAlpha','0.2','EdgeAlpha','0.2')
        text{count} = sprintf('%s p =< %.3f (proportion = %.3f)',method_type{nmethod},p_val_threshold(low_threshold(1)),mean(percent_shuffle_events{nmethod}{epoch}(low_threshold(1),:)));

        ylabel('Proportion of significant events')
        hold on
        xlabel('mean log odd difference')
        %             title(Behavioural_epoches{epoch})
        title('Original p value =< 0.05')
        count = count + 1;
    end
end
legend([s(1:4)], {text{1:4}})

xlim([-1 3])
ylim([0 0.7])

% Sig proportion and log odd at shuffle-corrected p value 0.05
subplot(2,2,2)
count = 1;
%         tempt = find(mean(percent_multi_events{nmethod}{nshuffle}{epoch},2) <= 0.05);
for nmethod = [7 8]
    for epoch = 2:3
        [c index(nmethod)] = min(abs(mean(percent_shuffle_events{nmethod}{epoch},2) - 0.05));
%          [c index(nmethod)] = min(abs(mean(percent_multi_events{nmethod}{1}{epoch},2) - 0.05));
        x = mean(mean(log_odd_difference{nmethod}{1}{epoch}(index(nmethod),:),2),2);
        y = mean(percent_sig_events{nmethod}{1}{epoch}(index(nmethod),:),2);
        s(count) = scatter(x,y,20,marker_shape{epoch},...
            'MarkerFaceColor',colour_line{nmethod},'MarkerEdgeColor',colour_line{nmethod},'MarkerFaceAlpha','0.5')
        %             hold on
        %             m{nshuffle} = errorbar(p_val_threshold(low_threshold),y,y-percent_multi_events_CI{nmethod}{nshuffle}{epoch}(low_threshold,1),...
        %                 percent_multi_events_CI{nmethod}{nshuffle}{epoch}(low_threshold,2)-y,'.','MarkerFaceColor',colour_line{nshuffle},'MarkerEdgeColor',colour_line{nshuffle})
        %             m{nshuffle}.Color = colour_line{nshuffle};
        hold on
        error_x = [log_odd_difference_CI{nmethod}{1}{epoch}(index(nmethod),1) log_odd_difference_CI{nmethod}{1}{epoch}(index(nmethod),2)...
            log_odd_difference_CI{nmethod}{1}{epoch}(index(nmethod),2) log_odd_difference_CI{nmethod}{1}{epoch}(index(nmethod),1)]
        error_y = [percent_sig_events_CI{nmethod}{1}{epoch}(index(nmethod),1) percent_sig_events_CI{nmethod}{1}{epoch}(index(nmethod),1)...
            percent_sig_events_CI{nmethod}{1}{epoch}(index(nmethod),2) percent_sig_events_CI{nmethod}{1}{epoch}(index(nmethod),2)]

        patch(error_x,error_y,colour_line{nmethod},'EdgeColor',colour_line{nmethod},'FaceAlpha','0.2','EdgeAlpha','0.2')
        text{count} = sprintf('%s p =< %.3f (proportion = %.3f)',method_type{nmethod},p_val_threshold(index(nmethod)),mean(percent_shuffle_events{nmethod}{epoch}(index(nmethod),:)));
%         text{nmethod} = sprintf('%s p =< %.3f (proportion = %.3f)',shuffle_type{nmethod},p_val_threshold(index(nmethod)),mean(percent_multi_events{nmethod}{1}{epoch}(index(nmethod),:)));

        ylabel('Proportion of significant events')
        hold on
        xlabel('mean log odd difference')
        %             title(Behavioural_epoches{epoch})
        title('Equivalent p value when cell-id shuffled sig proportion = 0.05')
        count = count + 1;
    end
end
%         plot([0 0.1],[0 0.3],'r')
legend([s(1:4)], {text{1:4}})
xlim([-1 3])
ylim([0 0.7])

cd ground_truth_original\Figure
filename = sprintf('linear optimisation original vs shuffle-corrected p value.pdf')
saveas(gcf,filename)
filename = sprintf('linear optimisation original vs shuffle-corrected p value.fig')
saveas(gcf,filename)
cd ..
cd ..
clf



%% All methods in one plot

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
randomisation_type = {'original','Cell ID randomisation','cross experiment','place field randomisation','spike train randomisation'};


% Inset for false positive rate at p value 0.05 and p value for false
% positive rate 0.05

for nshuffle = [2 4 5]

    fig = figure(nshuffle)
    fig.Position = [834 116 575 531];
    count = 1;

    subplot(2,2,1)
    for nmethod = 1:length(method_type)
        for epoch = 1:3
            mean_false_positive_rate(count) = mean(percent_shuffle_events{nmethod}{epoch}{nshuffle}(low_threshold(1),:),2);
            %         b = prctile(percent_shuffle_events_POST_ripple{nmethod}{epoch}(low_threshold(1),:),[2.5 97.5])

            s(nmethod) = scatter(epoch-0.5,mean_false_positive_rate(count),20,marker_shape{epoch},...
                'MarkerFaceColor',colour_line{nmethod},'MarkerEdgeColor',colour_line{nmethod},'MarkerFaceAlpha','0.5')
            hold on
            m{nmethod} = errorbar(epoch-0.5,mean_false_positive_rate(count),mean_false_positive_rate(count)-percent_shuffle_events_CI{nmethod}{epoch}{nshuffle}(low_threshold(1),1),...
                percent_shuffle_events_CI{nmethod}{epoch}{nshuffle}(low_threshold(1),2)-mean_false_positive_rate(count),'.','MarkerFaceColor',colour_line{nmethod},'MarkerEdgeColor',colour_line{nmethod})
            m{nmethod}.Color = colour_line{nmethod};
            %         set(gca,'yscale','log')
            count = count + 1;
        end
    end
    ylabel('mean false positive rate')
    pbaspect([1 3 1])
    ylim([0 0.2])
    xlim([0 3])
    set(gca,'xtick',[]);
    ax = gca;
    ax.FontSize = 12;
    set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);

    subplot(2,2,2)
    count = 1;
    for nmethod = 1:length(method_type)
        for epoch = 1:3
            [xx index] = min(abs(mean(percent_shuffle_events{nmethod}{epoch}{nshuffle},2) - 0.05));
            adjusted_pvalue(count) = p_val_threshold(index);

            s(nmethod) = scatter(epoch-0.5,adjusted_pvalue(count),20,marker_shape{epoch},...
                'MarkerFaceColor',colour_line{nmethod},'MarkerEdgeColor',colour_line{nmethod},'MarkerFaceAlpha','0.5')
            hold on
            %         m{nmethod} = errorbar(p_val_threshold(low_threshold),y,y-percent_shuffle_events_POST_ripple_CI{nmethod}{epoch}(low_threshold,1),...
            %             percent_shuffle_events_POST_ripple_CI{nmethod}{epoch}(low_threshold,2)-y,'.','MarkerFaceColor',colour_line{nmethod},'MarkerEdgeColor',colour_line{nmethod})
            %         m{nmethod}.Color = colour_line{nmethod};
            %         set(gca,'yscale','log')
            count = count + 1;
        end
    end
    ylabel('Adjusted p value')
    pbaspect([1 3 1])
    ylim([0 0.21])
    xlim([0 3])
    ax = gca;
    ax.FontSize = 12;
    set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
    set(gca,'xtick',[]);
    sgtitle(randomisation_type{nshuffle})

    cd ground_truth_original\Figure
    filename = sprintf('all methods in one false positive rate and p value inset (%s).pdf',randomisation_type{nshuffle})
    saveas(gcf,filename)
    filename = sprintf('all methods in one false positive rate and p value inset (%s).fig',randomisation_type{nshuffle})
    saveas(gcf,filename)
    cd ..
    cd ..
    clf
end


for nshuffle = [2 4 5]

    fig = figure(nshuffle)
    fig.Position = [834 116 850 885];
    subplot(2,2,1)
    count = 1;
    for nmethod = 1:length(method)
        for epoch = 1:3
            x = mean(mean(log_odd_difference{nmethod}{1}{epoch}(low_threshold(1),:),2),2);
            y = mean(percent_sig_events{nmethod}{1}{epoch}(low_threshold(1),:),2);
            s(count) = scatter(x,y,20,marker_shape{epoch},...
                'MarkerFaceColor',colour_line{nmethod},'MarkerEdgeColor',colour_line{nmethod},'MarkerFaceAlpha','0.5')
            %             hold on
            %             m{nshuffle} = errorbar(p_val_threshold(low_threshold),y,y-percent_multi_events_CI{nmethod}{nshuffle}{epoch}(low_threshold,1),...
            %                 percent_multi_events_CI{nmethod}{nshuffle}{epoch}(low_threshold,2)-y,'.','MarkerFaceColor',colour_line{nshuffle},'MarkerEdgeColor',colour_line{nshuffle})
            %             m{nshuffle}.Color = colour_line{nshuffle};
            hold on
            error_x = [log_odd_difference_CI{nmethod}{1}{epoch}(low_threshold(1),1) log_odd_difference_CI{nmethod}{1}{epoch}(low_threshold(1),2)...
                log_odd_difference_CI{nmethod}{1}{epoch}(low_threshold(1),2) log_odd_difference_CI{nmethod}{1}{epoch}(low_threshold(1),1)]
            error_y = [percent_sig_events_CI{nmethod}{1}{epoch}(low_threshold(1),1) percent_sig_events_CI{nmethod}{1}{epoch}(low_threshold(1),1)...
                percent_sig_events_CI{nmethod}{1}{epoch}(low_threshold(1),2) percent_sig_events_CI{nmethod}{1}{epoch}(low_threshold(1),2)]

            patch(error_x,error_y,colour_line{nmethod},'EdgeColor',colour_line{nmethod},'FaceAlpha','0.2','EdgeAlpha','0.2')
            text{count} = sprintf('%s p =< %.3f (proportion = %.3f)',method_type{nmethod},p_val_threshold(low_threshold(1)),mean(percent_shuffle_events{nmethod}{epoch}{nshuffle}(low_threshold(1),:)));

            ylabel('Proportion of significant events')
            hold on
            xlabel('mean log odd difference')
            %             title(Behavioural_epoches{epoch})
            title('Original p value =< 0.05')
            count = count + 1;
        end
        ax = gca;
        ax.FontSize = 12;
        set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
    end
    legend([s(1:24)],...
        {text{1:24}})

    xlim([-1 3])
    ylim([0 0.6])

    % Sig proportion and log odd at shuffle-corrected p value 0.05
    subplot(2,2,2)
    count = 1;
    %         tempt = find(mean(percent_multi_events{nmethod}{nshuffle}{epoch},2) <= 0.05);
    for nmethod = 1:length(method)
        for epoch = 1:3
            [c index(nmethod)] = min(abs(mean(percent_shuffle_events{nmethod}{epoch}{nshuffle},2) - 0.05));
            %          [c index(nmethod)] = min(abs(mean(percent_multi_events{nmethod}{1}{epoch},2) - 0.05));
            x = mean(mean(log_odd_difference{nmethod}{1}{epoch}(index(nmethod),:),2),2);
            y = mean(percent_sig_events{nmethod}{1}{epoch}(index(nmethod),:),2);
            s(count) = scatter(x,y,20,marker_shape{epoch},...
                'MarkerFaceColor',colour_line{nmethod},'MarkerEdgeColor',colour_line{nmethod},'MarkerFaceAlpha','0.5')
            %             hold on
            %             m{nshuffle} = errorbar(p_val_threshold(low_threshold),y,y-percent_multi_events_CI{nmethod}{nshuffle}{epoch}(low_threshold,1),...
            %                 percent_multi_events_CI{nmethod}{nshuffle}{epoch}(low_threshold,2)-y,'.','MarkerFaceColor',colour_line{nshuffle},'MarkerEdgeColor',colour_line{nshuffle})
            %             m{nshuffle}.Color = colour_line{nshuffle};
            hold on
            error_x = [log_odd_difference_CI{nmethod}{1}{epoch}(index(nmethod),1) log_odd_difference_CI{nmethod}{1}{epoch}(index(nmethod),2)...
                log_odd_difference_CI{nmethod}{1}{epoch}(index(nmethod),2) log_odd_difference_CI{nmethod}{1}{epoch}(index(nmethod),1)]
            error_y = [percent_sig_events_CI{nmethod}{1}{epoch}(index(nmethod),1) percent_sig_events_CI{nmethod}{1}{epoch}(index(nmethod),1)...
                percent_sig_events_CI{nmethod}{1}{epoch}(index(nmethod),2) percent_sig_events_CI{nmethod}{1}{epoch}(index(nmethod),2)]

            patch(error_x,error_y,colour_line{nmethod},'EdgeColor',colour_line{nmethod},'FaceAlpha','0.2','EdgeAlpha','0.2')
            text{count} = sprintf('%s p =< %.3f (proportion = %.3f)',method_type{nmethod},p_val_threshold(index(nmethod)),mean(percent_shuffle_events{nmethod}{epoch}{nshuffle}(index(nmethod),:)));
            %         text{nmethod} = sprintf('%s p =< %.3f (proportion = %.3f)',shuffle_type{nmethod},p_val_threshold(index(nmethod)),mean(percent_multi_events{nmethod}{1}{epoch}(index(nmethod),:)));

            ylabel('Proportion of significant events')
            hold on
            xlabel('mean log odd difference')
            %             title(Behavioural_epoches{epoch})
            title('Equivalent p value when cell-id shuffled sig proportion = 0.05')
            count = count + 1;
        end
        ax = gca;
        ax.FontSize = 12;
        set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
    end
    %         plot([0 0.1],[0 0.3],'r')
    legend([s(1:24)],...
        {text{1:24}})

    xlim([-1 3])
    ylim([0 0.6])
    sgtitle(randomisation_type{nshuffle})
end

cd ground_truth_original\Figure
filename = sprintf('all methods in one original vs shuffle-corrected p value.pdf')
saveas(gcf,filename)
filename = sprintf('all methods in one original vs shuffle-corrected p value.fig')
saveas(gcf,filename)
cd ..
cd ..
clf


%% All methods in one plot (shuffle-subtracted version)

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

for nshuffle = [2 4 5]

    fig = figure(nshuffle)
%     fig = figure(1)
    fig.Position = [834 116 850 885];
    subplot(2,2,1)
    count = 1;
    for nmethod = 1:length(method)
        for epoch = 1:3
            % shuffle-subtracted mean log odd difference
            shuffle_subtracted_log_odd_difference_CI{nmethod}{epoch} = prctile(log_odd_difference{nmethod}{1}{epoch} - log_odd_difference{nmethod}{2}{epoch},[2.5 97.5],2);
            shuffle_subtracted_log_odd_difference{nmethod}{epoch} = mean(log_odd_difference{nmethod}{1}{epoch} - log_odd_difference{nmethod}{2}{epoch},2);

            x_CI = shuffle_subtracted_log_odd_difference_CI{nmethod}{epoch}(low_threshold(1),:);
            x = shuffle_subtracted_log_odd_difference{nmethod}{epoch}(low_threshold(1),:);

            y = mean(percent_sig_events{nmethod}{1}{epoch}(low_threshold(1),:),2);
            s(count) = scatter(x,y,20,marker_shape{epoch},...
                'MarkerFaceColor',colour_line{nmethod},'MarkerEdgeColor',colour_line{nmethod},'MarkerFaceAlpha','0.5')
            %             hold on
            %             m{nshuffle} = errorbar(p_val_threshold(low_threshold),y,y-percent_multi_events_CI{nmethod}{nshuffle}{epoch}(low_threshold,1),...
            %                 percent_multi_events_CI{nmethod}{nshuffle}{epoch}(low_threshold,2)-y,'.','MarkerFaceColor',colour_line{nshuffle},'MarkerEdgeColor',colour_line{nshuffle})
            %             m{nshuffle}.Color = colour_line{nshuffle};
            hold on
            error_x = [x_CI(1) x_CI(2) x_CI(2) x_CI(1)];
            error_y = [percent_sig_events_CI{nmethod}{1}{epoch}(low_threshold(1),1) percent_sig_events_CI{nmethod}{1}{epoch}(low_threshold(1),1)...
                percent_sig_events_CI{nmethod}{1}{epoch}(low_threshold(1),2) percent_sig_events_CI{nmethod}{1}{epoch}(low_threshold(1),2)]

            patch(error_x,error_y,colour_line{nmethod},'EdgeColor',colour_line{nmethod},'FaceAlpha','0.2','EdgeAlpha','0.2')
            text{count} = sprintf('%s p =< %.3f (proportion = %.3f)',method_type{nmethod},p_val_threshold(low_threshold(1)),mean(percent_shuffle_events{nmethod}{epoch}{nshuffle}(low_threshold(1),:)));

            ylabel('Proportion of significant events')
            hold on
            xlabel('shuffle-subtracted mean log odd difference')
            %             title(Behavioural_epoches{epoch})
            title('Original p value =< 0.05')

            count = count + 1;
        end
        ax = gca;
        ax.FontSize = 12;
        set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
    end
    legend([s(1:24)],...
        {text{1:24}})

    xlim([-1 3])
    ylim([0 0.6])


    % Sig proportion and log odd at shuffle-corrected p value 0.05
    subplot(2,2,2)
    count = 1;
    %         tempt = find(mean(percent_multi_events{nmethod}{nshuffle}{epoch},2) <= 0.05);
    for nmethod = 1:length(method)
        for epoch = 1:3
            [c index(nmethod)] = min(abs(mean(percent_shuffle_events{nmethod}{epoch}{nshuffle},2) - 0.05));
            %          [c index(nmethod)] = min(abs(mean(percent_multi_events{nmethod}{1}{epoch},2) - 0.05));
            x_CI = shuffle_subtracted_log_odd_difference_CI{nmethod}{epoch}(index(nmethod),:);
            x = shuffle_subtracted_log_odd_difference{nmethod}{epoch}(index(nmethod),:);

            y = mean(percent_sig_events{nmethod}{1}{epoch}(index(nmethod),:),2);
            s(count) = scatter(x,y,20,marker_shape{epoch},...
                'MarkerFaceColor',colour_line{nmethod},'MarkerEdgeColor',colour_line{nmethod},'MarkerFaceAlpha','0.5')
            %             hold on
            %             m{nshuffle} = errorbar(p_val_threshold(low_threshold),y,y-percent_multi_events_CI{nmethod}{nshuffle}{epoch}(low_threshold,1),...
            %                 percent_multi_events_CI{nmethod}{nshuffle}{epoch}(low_threshold,2)-y,'.','MarkerFaceColor',colour_line{nshuffle},'MarkerEdgeColor',colour_line{nshuffle})
            %             m{nshuffle}.Color = colour_line{nshuffle};
            hold on
            error_x = [x_CI(1) x_CI(2) x_CI(2) x_CI(1)];
            error_y = [percent_sig_events_CI{nmethod}{1}{epoch}(index(nmethod),1) percent_sig_events_CI{nmethod}{1}{epoch}(index(nmethod),1)...
                percent_sig_events_CI{nmethod}{1}{epoch}(index(nmethod),2) percent_sig_events_CI{nmethod}{1}{epoch}(index(nmethod),2)]

            patch(error_x,error_y,colour_line{nmethod},'EdgeColor',colour_line{nmethod},'FaceAlpha','0.2','EdgeAlpha','0.2')
            text{count} = sprintf('%s p =< %.3f (proportion = %.3f)',method_type{nmethod},p_val_threshold(index(nmethod)),mean(percent_shuffle_events{nmethod}{epoch}{nshuffle}(index(nmethod),:)));
            %         text{nmethod} = sprintf('%s p =< %.3f (proportion = %.3f)',shuffle_type{nmethod},p_val_threshold(index(nmethod)),mean(percent_multi_events{nmethod}{1}{epoch}(index(nmethod),:)));

            ylabel('Proportion of significant events')
            hold on
            xlabel('shuffle-subtracted mean log odd difference')
            %             title(Behavioural_epoches{epoch})
            title('Equivalent p value when cell-id shuffled sig proportion = 0.05')
            count = count + 1;
        end
        ax = gca;
        ax.FontSize = 12;
        set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
    end
    %         plot([0 0.1],[0 0.3],'r')
    legend([s(1:24)],...
        {text{1:24}})

    xlim([-1 3])
    ylim([0 0.6])
    sgtitle(randomisation_type{nshuffle})
end

cd ground_truth_original\Figure
filename = sprintf('all methods in one (shuffle-corrected) original vs shuffle-corrected p value.pdf')
saveas(gcf,filename)
filename = sprintf('all methods in one (shuffle-corrected) original vs shuffle-corrected p value.fig')
saveas(gcf,filename)
cd ..
cd ..
clf

%% Replay score vs mean log odds difference
% p_val_threshold = round(10.^[-3:0.01:log10(0.2)],4);
% % p_val_threshold = 0.001:0.0001:0.2;
% p_val_threshold = flip(unique(p_val_threshold));
% p_val_threshold(1) = 0.2;
% low_threshold = find(ismembertol(p_val_threshold,[0.0501 0.02 0.01 0.005 0.002 0.001],eps) == 1);
% method_type = {'wcorr 1 shuffle','wcorr 1 shuffle + jump distance','wcorr 2 shuffles','wcorr 3 shuffles'...
%     'spearman median spike','spearman all spikes',...
%     'linear 1 shuffle','linear 2 shuffles'};
% 
% % 
% % colour_line= {[34,94,168]/255,...
% %     [253,141,60]/255,[254,217,11]/255,...
% %     [227,26,28]/255};
% colour_line= {[0,0,0]/255,[67,147,195]/255,[244,109,67]/255,[215,48,39]/255};
% alpha_level = linspace(0.2,1,length(low_threshold));
% for nmethod = 3:3
%     fig = figure(nmethod)
%     fig.Position = [834 116 850 700];
%     for epoch = 2:3
%         subplot(2,2,epoch-1)
%         for nshuffle = 1:2
% 
%             y = mean(replay_score_compare{nmethod}{nshuffle}{epoch},2)';
%             x = mean(log_odd_difference{nmethod}{nshuffle}{epoch},2)';
%             if nshuffle == 1
%                 for threshold = 1:length(low_threshold)
%                     scatter(x(low_threshold(threshold)),y(low_threshold(threshold)),'filled','MarkerFaceColor',colour_line{epoch},'MarkerFaceAlpha',alpha_level(threshold),'MarkerEdgeColor',colour_line{epoch})
%                     hold on
%                 end
%             elseif nshuffle == 2
%                 for threshold = 1:length(low_threshold)
%                     scatter(x(low_threshold(threshold)),y(low_threshold(threshold)),'filled','MarkerFaceColor','k','MarkerFaceAlpha',alpha_level(threshold),'MarkerEdgeColor',colour_line{epoch})
%                     hold on
%                 end
%             end
% 
%             UCI = log_odd_difference_CI{nmethod}{nshuffle}{epoch}(:,2)';
%             UCI(isnan(x)) = [];
%             LCI = log_odd_difference_CI{nmethod}{nshuffle}{epoch}(:,1)';
%             LCI(isnan(x)) = [];
%             y(isnan(x)) = [];
%             x(isnan(x)) = [];
%             if nshuffle == 1
%                 p(nshuffle) = patch([LCI fliplr(UCI)], [y fliplr(y)], colour_line{epoch},'FaceAlpha','0.2','LineStyle','none');
%                 hold on
%             elseif nshuffle == 2
%                 p(nshuffle) = patch([LCI fliplr(UCI)], [y fliplr(y)],'k','FaceAlpha','0.2','LineStyle','none');
%                 hold on
%             end
% 
%             title(Behavioural_epoches{epoch});
% 
%             xlabel('mean log odds difference')
%             ylabel('mean replay score')
%             ax = gca;
%             ax.FontSize = 12;
%             set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
%         end
%         %         legend([p(1),p(2),p(3)], {Behavioural_epoches{1},Behavioural_epoches{2},Behavioural_epoches{3}},'Position',[0.7 0.3 0.05 0.05])
% 
%     end
%     legend([p(1),p(2)], {'original data','shuffled data'},'Location','southeast')
%     %         legend([p(1),p(2)], {'original data','shuffled data'},'Position',[0.7 0.3 0.05 0.05])
%     cd ground_truth_original\Figure
%     filename = sprintf('%s replay score comparisions CI.pdf',method{nmethod})
%     sgtitle(method{nmethod})
%     saveas(gcf,filename)
%     cd ..
%     cd ..
%     clf
% end

colour_line= {[0,0,0]/255,[67,147,195]/255,[244,109,67]/255,[215,48,39]/255};
alpha_level = linspace(0.2,1,length(low_threshold));
for nmethod = 3:3
    fig = figure(nmethod)
    fig.Position = [834 116 850 700];
    for epoch = 2:3
        subplot(2,2,epoch-1)
        for nshuffle = 1:2

            y = mean(replay_score_compare{nmethod}{nshuffle}{epoch},2)';
            x = mean(log_odd_difference{nmethod}{nshuffle}{epoch},2)';
            if nshuffle == 1
                for threshold = 1:length(low_threshold)
                    scatter(x(low_threshold(threshold)),y(low_threshold(threshold)),'filled','MarkerFaceColor',colour_line{epoch},'MarkerFaceAlpha',alpha_level(threshold),'MarkerEdgeColor',colour_line{epoch})
                    hold on
                end
            elseif nshuffle == 2
                for threshold = 1:length(low_threshold)
                    scatter(x(low_threshold(threshold)),y(low_threshold(threshold)),'filled','MarkerFaceColor','k','MarkerFaceAlpha',alpha_level(threshold),'MarkerEdgeColor',colour_line{epoch})
                    hold on
                end
            end

            UCI = log_odd_difference_CI{nmethod}{nshuffle}{epoch}(:,2)';
            UCI(isnan(x)) = [];
            LCI = log_odd_difference_CI{nmethod}{nshuffle}{epoch}(:,1)';
            LCI(isnan(x)) = [];
            y(isnan(x)) = [];
            x(isnan(x)) = [];
            if nshuffle == 1
                p(nshuffle) = patch([LCI fliplr(UCI)], [y fliplr(y)], colour_line{epoch},'FaceAlpha','0.2','LineStyle','none');
                hold on
            elseif nshuffle == 2
                p(nshuffle) = patch([LCI fliplr(UCI)], [y fliplr(y)],'k','FaceAlpha','0.2','LineStyle','none');
                hold on
            end

            title(Behavioural_epoches{epoch});

            xlabel('mean log odds difference')
            ylabel('mean replay score')

        end
        ax = gca;
        ax.FontSize = 12;
        set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
        %         legend([p(1),p(2),p(3)], {Behavioural_epoches{1},Behavioural_epoches{2},Behavioural_epoches{3}},'Position',[0.7 0.3 0.05 0.05])
        legend([p(1),p(2)], {'original data','shuffled data'},'Location','southeast')
        %         legend([p(1),p(2)], {'original data','shuffled data'},'Position',[0.7 0.3 0.05 0.05])

    end
    cd ground_truth_original\Figure
    filename = sprintf('%s replay score comparisions CI.pdf',method{nmethod})
    sgtitle(method{nmethod})
    saveas(gcf,filename)
    cd ..
    cd ..
    clf
end


%% Cell id shuffle vs original 

p_val_threshold = round(10.^[-3:0.01:log10(0.2)],4);
% p_val_threshold = 0.001:0.0001:0.2;
p_val_threshold = flip(unique(p_val_threshold));
p_val_threshold(1) = 0.2;
low_threshold = find(ismembertol(p_val_threshold,[0.0501 0.02 0.01 0.005 0.002 0.001],eps) == 1);


method_type = {'wcorr 1 shuffle','wcorr 1 shuffle + jump distance','wcorr 2 shuffles','wcorr 3 shuffles'...
    'spearman median spike','spearman all spikes',...
    'linear 1 shuffle','linear 2 shuffles'};

colour_line= {[127,205,187]/255,[65,182,196]/255,[34,94,168]/255,[37,52,148]/255,...
    [253,141,60]/255,[254,217,11]/255,...
    [128,0,38]/255,[227,26,28]/255};

alpha_level = linspace(0.2,0.7,length(low_threshold));

for nmethod = 1:length(method)
    for epoch = 1:3
        fig = figure(2)
        fig.Position = [834 116 850 700];
        subplot(2,2,epoch)

        for nshuffle = 1:2
            x = mean(log_odd_difference{nmethod}{nshuffle}{epoch},2)';
            y = mean(percent_sig_events{nmethod}{nshuffle}{epoch},2)';

            if nshuffle == 1
                for threshold = 1:length(low_threshold)
                    scatter(x(low_threshold(threshold)),y(low_threshold(threshold)),'filled','MarkerFaceColor',colour_line{nmethod},'MarkerFaceAlpha',alpha_level(threshold),'MarkerEdgeColor',colour_line{nmethod})
                    hold on
                end
            elseif nshuffle == 2
                for threshold = 1:length(low_threshold)
                    scatter(x(low_threshold(threshold)),y(low_threshold(threshold)),'filled','MarkerFaceColor','k','MarkerFaceAlpha',alpha_level(threshold),'MarkerEdgeColor',colour_line{nmethod})
                    hold on
                end
            end

            UCI = log_odd_difference_CI{nmethod}{nshuffle}{epoch}(:,2)';
            UCI(isnan(x)) = [];
            LCI = log_odd_difference_CI{nmethod}{nshuffle}{epoch}(:,1)';
            LCI(isnan(x)) = [];
            y(isnan(x)) = [];
            x(isnan(x)) = [];
            if nshuffle == 1
                p(nshuffle) = patch([LCI fliplr(UCI)], [y fliplr(y)], colour_line{nmethod},'FaceAlpha','0.2','LineStyle','none');
                hold on
            elseif nshuffle == 2
                p(nshuffle) = patch([LCI fliplr(UCI)], [y fliplr(y)],'k','FaceAlpha','0.2','LineStyle','none');
                hold on
            end
            xlabel('mean log odds difference')
            ylabel('Proportion of all replay events')
            title(Behavioural_epoches{epoch});
        end
    end
%     ylim([0 1])
    ax = gca;
    ax.FontSize = 12;
    set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
    legend([p(1),p(2)], {'Original','Cell-id shuffle'},'Position',[0.7 0.3 0.05 0.05])
    
    cd ground_truth_original\Figure
    filename = sprintf('%s optimisation comparisions CI.pdf',method{nmethod})
    sgtitle(method{nmethod})
    saveas(gcf,filename)
    cd ..
    cd ..
    clf
end
%     legend([p(1) p(3) p(6)], {'Shuffle 3','Shuffle 2+3','Shuffle 1+2+3'})





%% Method optimisation
%     shuffle_type = {'POST Place','POST place + POST time','PRE place + POST time','PRE place + POST place','PRE place + PRE spike','PRE place + PRE spike + POST Place'}';
% ripple_type = {'Ripple Threshold 1','Ripple Threshold 3','Ripple Threshold 5','Ripple Threshold 10'}';

method_type = {'wcorr 1 shuffle','wcorr 1 shuffle + jump distance','wcorr 2 shuffles','wcorr 3 shuffles'...
    'spearman median spike','spearman all spikes',...
    'linear 1 shuffle','linear 2 shuffles'};

low_threshold = find(ismembertol(p_val_threshold,[0.0501 0.02 0.01 0.005 0.002 0.001],eps) == 1);
shuffle_name(1) = {'Shuffle Type'};
behave_state(1) = {'Behavioural State'};

Mean_5(1) = {'Mean Log Odd difference at p =< 0.05'};
Mean_e(1) = {'Mean Log Odd difference at equivalent p value'};

UCI_5(1) = {'Upper CI at p =< 0.05'};
UCI_e(1) = {'Upper CI at equivalent p value'};

LCI_5(1) = {'Lower CI at p =< 0.05'};
LCI_e(1) = {'Lower CI at equivalent p value'};

proportion_5(1) = {'Proportion at p =< 0.05'};
proportion_5s(1) = {'Shuffle propotion at p =< 0.05'};
proportion_e(1) = {'Equivalent significant propotion'};
proportion_es(1) = {'Equivalent shuffle propotion'};

LCI_proportion_5(1) = {'Lower CI Proportion at p =< 0.05'};
UCI_proportion_5(1) = {'Upper CI Proportion at p =< 0.05'};
LCI_proportion_5s(1) = {'Lower CI shuffle Proportion at p =< 0.05'};
UCI_proportion_5s(1) = {'Upper CI shuffle propotion at p =< 0.05'};
LCI_proportion_e(1) = {'Lower CI equivalent significant proportion'};
UCI_proportion_e(1) = {'Upper CI equivalent significant proportion'};
LCI_proportion_es(1) = {'Lower CI Equivalent shuffle propotion'};
UCI_proportion_es(1) = {'Upper CI Equivalent shuffle propotion'};

pvalue_e(1) = {'p value at equivalent significant proportion'};

c = 2;

for epoch = 1:3
    %     tempt = [mean(percent_sig_events{nmethod}{1}{epoch}(low_threshold(1),:)) mean(percent_sig_events{nmethod}{2}{epoch}(low_threshold(1),:))...
    %         mean(percent_sig_events{nmethod}{3}{epoch}(low_threshold(1),:)) mean(percent_sig_events{nmethod}{4}{epoch}(low_threshold(1),:))];
    %     [minimum mindex] = min(tempt);
    for nmethod = 1:length(method)
        for nshuffle = 1:1
            shuffle_name(c) = {method_type{nmethod}};
            Mean_5(c) = {mean(log_odd_difference{nmethod}{nshuffle}{epoch}(low_threshold(1),:))};
            a = prctile(log_odd_difference{nmethod}{nshuffle}{epoch}(low_threshold(1),:),[2.5 97.5]);
            LCI_5(c) = {a(1)};
            UCI_5(c) = {a(2)};

            proportion_5(c) = {mean(percent_sig_events{nmethod}{nshuffle}{epoch}(low_threshold(1),:))};
            b = prctile(percent_sig_events{nmethod}{nshuffle}{epoch}(low_threshold(1),:),[2.5 97.5]);
            LCI_proportion_5(c) = {b(1)};
            UCI_proportion_5(c) = {b(2)};

            proportion_5s(c) = {mean(percent_shuffle_events{nmethod}{epoch}(low_threshold(1),:))};
            b = prctile(percent_shuffle_events{nmethod}{epoch}(low_threshold(1),:),[2.5 97.5]);
            LCI_proportion_5s(c) = {b(1)};
            UCI_proportion_5s(c) = {b(2)};

            % Get p value index at shuffle-corrected equivalent propotion of events
            %         [xx index] = min(abs(mean(percent_sig_events{nmethod}{nshuffle}{epoch},2) - mean(percent_sig_events{nmethod}{mindex}{epoch}(low_threshold(1),:))));
            [xx index] = min(abs(mean(percent_shuffle_events{nmethod}{epoch},2) - 0.05));
            Mean_e(c) = {mean(log_odd_difference{nmethod}{nshuffle}{epoch}(index,:))};
            a = prctile(log_odd_difference{nmethod}{nshuffle}{epoch}(index,:),[2.5 97.5]);
            LCI_e(c) = {a(1)};
            UCI_e(c) = {a(2)};
            %         proportion_e(c) = {mean(percent_sig_events{nmethod}{nshuffle}{epoch}(index,:))};
            proportion_e(c) = {mean(percent_sig_events{nmethod}{nshuffle}{epoch}(index,:))};
            %         b = prctile(percent_sig_events{nmethod}{nshuffle}{epoch}(index,:),[2.5 97.5]);
            b = prctile(percent_sig_events{nmethod}{nshuffle}{epoch}(index,:),[2.5 97.5]);
            LCI_proportion_e(c) = {b(1)};
            UCI_proportion_e(c) = {b(2)};

            proportion_es(c) = {mean(percent_shuffle_events{nmethod}{epoch}(index,:))};
            b = prctile(percent_shuffle_events{nmethod}{epoch}(index,:),[2.5 97.5]);
            LCI_proportion_es(c) = {b(1)};
            UCI_proportion_es(c) = {b(2)};
            pvalue_e(c) = {p_val_threshold(index)};

            if epoch == 1
                behave_state(c) = {'PRE'};
            elseif epoch == 2

                behave_state(c) = {'RUN'};
            elseif epoch == 3

                behave_state(c) = {'POST'};
            end
            c = c+1;
        end


    end
end

cd ground_truth_original\Figure\
Table = table(behave_state',shuffle_name',Mean_5',UCI_5',LCI_5',proportion_5',UCI_proportion_5',LCI_proportion_5',...
    proportion_5s',UCI_proportion_5s',LCI_proportion_5s',...
    Mean_e',UCI_e',LCI_e',proportion_e',UCI_proportion_e',LCI_proportion_5',...
    proportion_es',UCI_proportion_es',LCI_proportion_es',pvalue_e');

writetable(Table,'method optimisation table.xlsx')
cd ..
cd ..



%% Method comparision table (shuffle subtracted log odd)
%     shuffle_type = {'POST Place','POST place + POST time','PRE place + POST time','PRE place + POST place','PRE place + PRE spike','PRE place + PRE spike + POST Place'}';
% ripple_type = {'Ripple Threshold 1','Ripple Threshold 3','Ripple Threshold 5','Ripple Threshold 10'}';

method_type = {'wcorr 1 shuffle','wcorr 1 shuffle + jump distance','wcorr 2 shuffles','wcorr 3 shuffles'...
    'spearman median spike','spearman all spikes',...
    'linear 1 shuffle','linear 2 shuffles'};

low_threshold = find(ismembertol(p_val_threshold,[0.0501 0.02 0.01 0.005 0.002 0.001],eps) == 1);
shuffle_name(1) = {'Shuffle Type'};
behave_state(1) = {'Behavioural State'};

Mean_5(1) = {'Shuffle-subtracted mean log odd difference at p =< 0.05'};
Mean_e(1) = {'Shuffle-subtracted mean log Odd difference at equivalent p value'};

UCI_5(1) = {'Upper CI at p =< 0.05'};
UCI_e(1) = {'Upper CI at equivalent p value'};

LCI_5(1) = {'Lower CI at p =< 0.05'};
LCI_e(1) = {'Lower CI at equivalent p value'};

proportion_5(1) = {'Proportion at p =< 0.05'};
proportion_5s(1) = {'Shuffle propotion at p =< 0.05'};
proportion_e(1) = {'Equivalent significant propotion'};
proportion_es(1) = {'Equivalent shuffle propotion'};

LCI_proportion_5(1) = {'Lower CI Proportion at p =< 0.05'};
UCI_proportion_5(1) = {'Upper CI Proportion at p =< 0.05'};
LCI_proportion_5s(1) = {'Lower CI shuffle Proportion at p =< 0.05'};
UCI_proportion_5s(1) = {'Upper CI shuffle propotion at p =< 0.05'};
LCI_proportion_e(1) = {'Lower CI equivalent significant proportion'};
UCI_proportion_e(1) = {'Upper CI equivalent significant proportion'};
LCI_proportion_es(1) = {'Lower CI Equivalent shuffle propotion'};
UCI_proportion_es(1) = {'Upper CI Equivalent shuffle propotion'};

pvalue_e(1) = {'p value at equivalent significant proportion'};

c = 2;

for epoch = 1:3
    %     tempt = [mean(percent_sig_events{nmethod}{1}{epoch}(low_threshold(1),:)) mean(percent_sig_events{nmethod}{2}{epoch}(low_threshold(1),:))...
    %         mean(percent_sig_events{nmethod}{3}{epoch}(low_threshold(1),:)) mean(percent_sig_events{nmethod}{4}{epoch}(low_threshold(1),:))];
    %     [minimum mindex] = min(tempt);
    for nmethod = 1:length(method)
        for nshuffle = 1:1
            shuffle_name(c) = {method_type{nmethod}};
            Mean_5(c) = {shuffle_subtracted_log_odd_difference{nmethod}{epoch}(low_threshold(1))};
            LCI_5(c) = {shuffle_subtracted_log_odd_difference_CI{nmethod}{epoch}(low_threshold(1),1)};
            UCI_5(c) = {shuffle_subtracted_log_odd_difference_CI{nmethod}{epoch}(low_threshold(1),2)};

            proportion_5(c) = {mean(percent_sig_events{nmethod}{nshuffle}{epoch}(low_threshold(1),:))};
            b = prctile(percent_sig_events{nmethod}{nshuffle}{epoch}(low_threshold(1),:),[2.5 97.5]);
            LCI_proportion_5(c) = {b(1)};
            UCI_proportion_5(c) = {b(2)};

            proportion_5s(c) = {mean(percent_shuffle_events{nmethod}{epoch}(low_threshold(1),:))};
            b = prctile(percent_shuffle_events{nmethod}{epoch}(low_threshold(1),:),[2.5 97.5]);
            LCI_proportion_5s(c) = {b(1)};
            UCI_proportion_5s(c) = {b(2)};

            % Get p value index at shuffle-corrected equivalent propotion of events
            %         [xx index] = min(abs(mean(percent_sig_events{nmethod}{nshuffle}{epoch},2) - mean(percent_sig_events{nmethod}{mindex}{epoch}(low_threshold(1),:))));
            [xx index] = min(abs(mean(percent_shuffle_events{nmethod}{epoch},2) - 0.05));
            Mean_e(c) = {shuffle_subtracted_log_odd_difference{nmethod}{epoch}(index)};
            LCI_e(c) = {shuffle_subtracted_log_odd_difference_CI{nmethod}{epoch}(index,1)};
            UCI_e(c) = {shuffle_subtracted_log_odd_difference_CI{nmethod}{epoch}(index,2)};

            %         proportion_e(c) = {mean(percent_sig_events{nmethod}{nshuffle}{epoch}(index,:))};
            proportion_e(c) = {mean(percent_sig_events{nmethod}{nshuffle}{epoch}(index,:))};
            %         b = prctile(percent_sig_events{nmethod}{nshuffle}{epoch}(index,:),[2.5 97.5]);
            b = prctile(percent_sig_events{nmethod}{nshuffle}{epoch}(index,:),[2.5 97.5]);
            LCI_proportion_e(c) = {b(1)};
            UCI_proportion_e(c) = {b(2)};

            proportion_es(c) = {mean(percent_shuffle_events{nmethod}{epoch}(index,:))};
            b = prctile(percent_shuffle_events{nmethod}{epoch}(index,:),[2.5 97.5]);
            LCI_proportion_es(c) = {b(1)};
            UCI_proportion_es(c) = {b(2)};
            pvalue_e(c) = {p_val_threshold(index)};

            if epoch == 1
                behave_state(c) = {'PRE'};
            elseif epoch == 2

                behave_state(c) = {'RUN'};
            elseif epoch == 3

                behave_state(c) = {'POST'};
            end
            c = c+1;
        end


    end
end

cd ground_truth_original\Figure\
Table = table(behave_state',shuffle_name',Mean_5',UCI_5',LCI_5',proportion_5',UCI_proportion_5',LCI_proportion_5',...
    proportion_5s',UCI_proportion_5s',LCI_proportion_5s',...
    Mean_e',UCI_e',LCI_e',proportion_e',UCI_proportion_e',LCI_proportion_e',...
    proportion_es',UCI_proportion_es',LCI_proportion_es',pvalue_e');

writetable(Table,'shuffle-subtracted method comparision table.xlsx')
cd ..
cd ..


end



