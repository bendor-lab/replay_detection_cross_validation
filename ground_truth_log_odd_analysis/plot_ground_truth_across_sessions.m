function [] = plot_ground_truth_across_sessions(folders,option,method)

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
session_total_number{3} = session_total_number{2};
total_number{3} = total_number{2};
total_number{4} = total_number{2};


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
%         load log_odd_wcorr_place_bin_circular_shift_place_field_shifted
%         log_odd_compare{nmethod}{4} = log_odd;

    elseif strcmp(method{nmethod},'wcorr 1 shuffle + jump distance')
        load log_odd_wcorr_place_bin_circular_shift
        log_odd_compare{nmethod}{1} = log_odd;
        load log_odd_wcorr_place_bin_circular_shift_global_remapped
        log_odd_compare{nmethod}{2} = log_odd;
        load log_odd_wcorr_place_bin_circular_shift_cross_experiment_shuffled
        log_odd_compare{nmethod}{3} = log_odd;
%         load log_odd_wcorr_place_bin_circular_shift_place_field_shifted
%         log_odd_compare{nmethod}{4} = log_odd;

    elseif strcmp(method{nmethod},'wcorr 2 shuffles')
        load log_odd_wcorr_PRE_place_POST_time
        log_odd_compare{nmethod}{1} = log_odd;
        load log_odd_wcorr_PRE_place_POST_time_global_remapped
        log_odd_compare{nmethod}{2} = log_odd;
        load log_odd_wcorr_PRE_place_POST_time_cross_experiment_shuffled
        log_odd_compare{nmethod}{3} = log_odd;
%         load log_odd_wcorr_PRE_place_POST_time_place_field_shifted
%         log_odd_compare{nmethod}{4} = log_odd;

    elseif strcmp(method{nmethod},'wcorr 3 shuffles')
        load log_odd_wcorr
        log_odd_compare{nmethod}{1} = log_odd;
        load log_odd_wcorr_global_remapped
        log_odd_compare{nmethod}{2} = log_odd;
        load log_odd_wcorr_cross_experiment_shuffled
        log_odd_compare{nmethod}{3} = log_odd;
%         load log_odd_wcorr_place_field_shifted
%         log_odd_compare{nmethod}{4} = log_odd;

    elseif strcmp(method{nmethod},'spearman median spike')
        load log_odd_spearman
        log_odd_compare{nmethod}{1} = log_odd;
        load log_odd_spearman_global_remapped
        log_odd_compare{nmethod}{2} = log_odd;
        load log_odd_spearman_cross_experiment_shuffled
        log_odd_compare{nmethod}{3} = log_odd;
%         load log_odd_spearman_place_field_shifted
%         log_odd_compare{nmethod}{4} = log_odd;

    elseif strcmp(method{nmethod},'spearman all spikes')
        load log_odd_spearman_all_spikes
        log_odd_compare{nmethod}{1} = log_odd;
        load log_odd_spearman_all_spikes_global_remapped
        log_odd_compare{nmethod}{2} = log_odd;
        load log_odd_spearman_all_spikes_cross_experiment_shuffled
        log_odd_compare{nmethod}{3} = log_odd;
%         load log_odd_spearman_all_spikes_place_field_shifted
%         log_odd_compare{nmethod}{4} = log_odd;

    elseif strcmp(method{nmethod},'linear 1 shuffle')
        load log_odd_linear_place_bin_circular_shift
        log_odd_compare{nmethod}{1} = log_odd;
        load log_odd_linear_place_bin_circular_shift_global_remapped
        log_odd_compare{nmethod}{2} = log_odd;
        load log_odd_linear_place_bin_circular_shift_cross_experiment_shuffled
        log_odd_compare{nmethod}{3} = log_odd;
%         load log_odd_linear_place_bin_circular_shift_place_field_shifted
%         log_odd_compare{nmethod}{4} = log_odd;
        
    elseif strcmp(method{nmethod},'linear 2 shuffles')
        load log_odd_linear_PRE_place_POST_time
        log_odd_compare{nmethod}{1} = log_odd;
        load log_odd_linear_PRE_place_POST_time_global_remapped
        log_odd_compare{nmethod}{2} = log_odd;
        load log_odd_linear_PRE_place_POST_time_cross_experiment_shuffled
        log_odd_compare{nmethod}{3} = log_odd;
%         load log_odd_linear_PRE_place_POST_time_place_field_shifted
%         log_odd_compare{nmethod}{4} = log_odd;
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
%         data{nmethod}{4} = log_odd_compare{nmethod}{4}.common_zscored.place_field_shifted_original;


        for nshuffle = 1:length(log_odd_compare{nmethod})
            log_pval{nmethod}{nshuffle} = log_odd_compare{nmethod}{nshuffle}.pvalue;
            segment_id{nmethod}{nshuffle} = log_odd_compare{nmethod}{nshuffle}.segment_id;
            max_jump{nmethod}{nshuffle} = log_odd_compare{nmethod}{nshuffle}.max_jump;
            replay_score{nmethod}{nshuffle} = log_odd_compare{nmethod}{nshuffle}.best_score;
        end

    elseif strcmp(option,'original')
        data{nmethod}{1} = log_odd_compare{nmethod}{1}.normal_zscored.original;
        data{nmethod}{2} = log_odd_compare{nmethod}{2}.normal_zscored.global_remapped_original;
        data{nmethod}{3} = log_odd_compare{nmethod}{3}.normal_zscored.cross_experiment_shuffled_original;
%         data{nmethod}{4} = log_odd_compare{nmethod}{4}.normal_zscored.place_field_shifted_original;


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
if exist('log_odd_difference_across_sessions.mat', 'file') ~= 2
    for nsession = 1:10
        for nmethod = 1:length(method)
            for epoch = 1:3
                for nshuffle = 1:length(log_pval{nmethod})
                    tic

                    parfor threshold = 1:length(p_val_threshold)
                        for nboot = 1:1000 % Bootstrapping 1000 times
                            s = RandStream('mrg32k3a','Seed',nboot); % Set random seed for resampling

                            % Resample canditate events with replacement

                            session_events = intersect(epoch_index{nmethod}{nshuffle}{epoch},find(log_odd_compare{nmethod}{nshuffle}.experiment == nsession));
                            resampled_event = datasample(s,session_events,length(session_events));

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
                            boot_percent_sig_events(threshold,nboot) = (length(track_1_index)+length(track_2_index)-multi_event_number)/ session_total_number{nshuffle}(nsession,epoch);

                            if nshuffle ~= 1
                                % cell id shuffled mean significant event proportion
                                boot_percent_shuffle_events(threshold,nboot) = ((length(track_1_index)+length(track_2_index))/2) / session_total_number{nshuffle}(nsession,epoch);
                            end
                        end
                    end

                    log_odd_difference{nsession}{nmethod}{nshuffle}{epoch} = boot_log_odd_difference;
                    log_odd_difference_CI{nsession}{nmethod}{nshuffle}{epoch} = prctile(boot_log_odd_difference,[2.5 97.5],2);
                    percent_sig_events{nsession}{nmethod}{nshuffle}{epoch} = boot_percent_sig_events;
                    percent_sig_events_CI{nsession}{nmethod}{nshuffle}{epoch} = prctile(boot_percent_sig_events,[2.5 97.5],2);
                    percent_multi_events{nsession}{nmethod}{nshuffle}{epoch} = boot_percent_multi_events;
                    percent_multi_events_CI{nsession}{nmethod}{nshuffle}{epoch} = prctile(boot_percent_multi_events,[2.5 97.5],2);
                    replay_score_compare{nsession}{nmethod}{nshuffle}{epoch} = boot_replay_score;
                    replay_score_compare_CI{nsession}{nmethod}{nshuffle}{epoch} = prctile(boot_replay_score,[2.5 97.5],2);

                    if nshuffle ~= 1
                        percent_shuffle_events{nsession}{nmethod}{epoch}{nshuffle} = boot_percent_shuffle_events;
                        percent_shuffle_events_CI{nsession}{nmethod}{epoch}{nshuffle} = prctile(boot_percent_shuffle_events,[2.5 97.5],2);
                    end

                    toc
                end
            end

        end
    end
    
    save log_odd_difference_across_sessions log_odd_difference log_odd_difference_CI percent_sig_events percent_sig_events_CI...
        percent_multi_events percent_multi_events_CI percent_shuffle_events percent_shuffle_events_CI replay_score_compare replay_score_compare_CI
    clear boot_log_odd_difference boot_percent_sig_events boot_percent_multi_events boot_percent_shuffle_events boot_replay_score
    cd ..
else
    load log_odd_difference_across_sessions
    cd ..
end



%% across session comparision

p_val_threshold = round(10.^[-3:0.01:log10(0.2)],4);
p_val_threshold = flip(unique(p_val_threshold));
p_val_threshold(1) = 0.2;
low_threshold = find(ismembertol(p_val_threshold,[0.0501 0.02 0.01 0.005 0.002 0.001],eps) == 1);
Behavioural_epoches = {'PRE','RUN','POST'};
marker_shape = {'o','+','*','.','x','square','diamond','^','v','pentagram'};
nfig = 1;
method_type = {'wcorr 2 shuffles',...
    'spearman median spike','spearman all spikes'};
colour_line= {[0 0 0],[34,94,168]/255,[253,141,60]/255}
% colour_line= {[34,94,168]/255,...
%     [253,141,60]/255,[254,217,11]/255};
% Sig proportion and log odd at p value 0.05

for nmethod = 1:2
    count = 1;
    fig = figure(1)
    fig.Position = [834 116 850 885];

    for nsession = 1:10

        subplot(2,2,nmethod)

        for epoch = 1:3
            x = mean(mean(log_odd_difference{nsession}{nmethod}{1}{epoch}(low_threshold(1),:),2),2);
            y = mean(percent_sig_events{nsession}{nmethod}{1}{epoch}(low_threshold(1),:),2);
            s(count) = scatter(x,y,20,marker_shape{nsession},...
                'MarkerFaceColor',colour_line{epoch},'MarkerEdgeColor',colour_line{epoch},'MarkerFaceAlpha','0.5')
            %             hold on
            %             m{nshuffle} = errorbar(p_val_threshold(low_threshold),y,y-percent_multi_events_CI{nmethod}{nshuffle}{epoch}(low_threshold,1),...
            %                 percent_multi_events_CI{nmethod}{nshuffle}{epoch}(low_threshold,2)-y,'.','MarkerFaceColor',colour_line{nshuffle},'MarkerEdgeColor',colour_line{nshuffle})
            %             m{nshuffle}.Color = colour_line{nshuffle};
            hold on
            error_x = [log_odd_difference_CI{nsession}{nmethod}{1}{epoch}(low_threshold(1),1) log_odd_difference_CI{nsession}{nmethod}{1}{epoch}(low_threshold(1),2)...
                log_odd_difference_CI{nsession}{nmethod}{1}{epoch}(low_threshold(1),2) log_odd_difference_CI{nsession}{nmethod}{1}{epoch}(low_threshold(1),1)]
            error_y = [percent_sig_events_CI{nsession}{nmethod}{1}{epoch}(low_threshold(1),1) percent_sig_events_CI{nsession}{nmethod}{1}{epoch}(low_threshold(1),1)...
                percent_sig_events_CI{nsession}{nmethod}{1}{epoch}(low_threshold(1),2) percent_sig_events_CI{nsession}{nmethod}{1}{epoch}(low_threshold(1),2)]

            patch(error_x,error_y,colour_line{epoch},'EdgeColor',colour_line{epoch},'FaceAlpha','0.2','EdgeAlpha','0.2')
            text{count} = sprintf('%s p =< %.3f (proportion = %.3f)',method_type{nmethod},p_val_threshold(low_threshold(1)),mean(percent_shuffle_events{nsession}{nmethod}{epoch}{nshuffle}(low_threshold(1),:)));

            ylabel('Proportion of significant events')
            hold on
            xlabel('mean log odds difference')
            %             title(Behavioural_epoches{epoch})
            title('Original alpha level =< 0.05')

            count = count + 1;
        end
        ax = gca;
        ax.FontSize = 12;
        set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
    end
    legend([s(1:4)], {text{1:4}})

    % xlim([-1 3])
    ylim([0 0.7])
end
%         plot([0 0.1],[0 0.3],'r')

cd ground_truth_original\resubmission_figure
filename = sprintf('across sessions comparisions CI.pdf')
saveas(gcf,filename)
filename = sprintf('across sessions comparisions CI.fig')
saveas(gcf,filename)
cd ..
cd ..
clf


for nmethod = 1:3
    % Sig proportion and log odd at shuffle-corrected p value 0.05
    count = 1;
    %         tempt = find(mean(percent_multi_events{nmethod}{nshuffle}{epoch},2) <= 0.05);
    fig = figure(2)
    fig.Position = [834 116 850 885];
    for nsession = 1:10

        subplot(2,2,nmethod)

        for epoch = 1:3
            [c index(nmethod)] = min(abs(mean(percent_shuffle_events{nsession}{nmethod}{epoch}{nshuffle},2) - 0.05));
            %          [c index(nmethod)] = min(abs(mean(percent_multi_events{nmethod}{1}{epoch},2) - 0.05));
            x = mean(mean(log_odd_difference{nsession}{nmethod}{1}{epoch}(index(nmethod),:),2),2);
            y = mean(percent_sig_events{nsession}{nmethod}{1}{epoch}(index(nmethod),:),2);
            s(count) = scatter(x,y,20,marker_shape{nsession},...
                'MarkerFaceColor',colour_line{epoch},'MarkerEdgeColor',colour_line{epoch},'MarkerFaceAlpha','0.5')
            %             hold on
            %             m{nshuffle} = errorbar(p_val_threshold(low_threshold),y,y-percent_multi_events_CI{nmethod}{nshuffle}{epoch}(low_threshold,1),...
            %                 percent_multi_events_CI{nmethod}{nshuffle}{epoch}(low_threshold,2)-y,'.','MarkerFaceColor',colour_line{nshuffle},'MarkerEdgeColor',colour_line{nshuffle})
            %             m{nshuffle}.Color = colour_line{nshuffle};
            hold on
            error_x = [log_odd_difference_CI{nsession}{nmethod}{1}{epoch}(index(nmethod),1) log_odd_difference_CI{nsession}{nmethod}{1}{epoch}(index(nmethod),2)...
                log_odd_difference_CI{nsession}{nmethod}{1}{epoch}(index(nmethod),2) log_odd_difference_CI{nsession}{nmethod}{1}{epoch}(index(nmethod),1)]
            error_y = [percent_sig_events_CI{nsession}{nmethod}{1}{epoch}(index(nmethod),1) percent_sig_events_CI{nsession}{nmethod}{1}{epoch}(index(nmethod),1)...
                percent_sig_events_CI{nsession}{nmethod}{1}{epoch}(index(nmethod),2) percent_sig_events_CI{nsession}{nmethod}{1}{epoch}(index(nmethod),2)]

            patch(error_x,error_y,colour_line{epoch},'EdgeColor',colour_line{epoch},'FaceAlpha','0.2','EdgeAlpha','0.2')
            text{count} = sprintf('%s p =< %.3f (proportion = %.3f)',method_type{nmethod},p_val_threshold(index(nmethod)),mean(percent_shuffle_events{nsession}{nmethod}{epoch}{nshuffle}(index(nmethod),:)));
            %         text{nmethod} = sprintf('%s p =< %.3f (proportion = %.3f)',shuffle_type{nmethod},p_val_threshold(index(nmethod)),mean(percent_multi_events{nmethod}{1}{epoch}(index(nmethod),:)));

            ylabel('Proportion of significant events')
            hold on
            xlabel('mean log odds difference')
            %             title(Behavioural_epoches{epoch})
            title('Equivalent alpha level when mean false-positive rate = 0.05')
            count = count+ 1;
        end
    end
    ax = gca;
    ax.FontSize = 12;
    set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
    legend([s(1:4)], {text{1:4}})
    ylim([0 0.7])
end
% xlim([-1 3])
cd ground_truth_original\resubmission_figure
filename = sprintf('across sessions comparisions FPR-macthed CI.pdf')
saveas(gcf,filename)
filename = sprintf('across sessions comparisions FPR-macthed CI.fig')
saveas(gcf,filename)
cd ..
cd ..
clf



end



