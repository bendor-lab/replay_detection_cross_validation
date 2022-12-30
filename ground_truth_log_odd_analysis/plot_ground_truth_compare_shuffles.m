function [] = plot_ground_truth_compare_shuffles(folders,option,method)
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

            elseif event_time >= time_range.track(1).behaviour(1) & event_time <= time_range.track(1).behaviour(2) %If RUN Track 1
                total_number{nshuffle}(2) = total_number{nshuffle}(2) + 1*number_of_global_remapped_shuffles;
                session_total_number{nshuffle}(nfolder,2) = session_total_number{nshuffle}(nfolder,2) + 1*number_of_global_remapped_shuffles;

            elseif event_time >= time_range.track(2).behaviour(1) & event_time <= time_range.track(2).behaviour(2) %If RUN Track 2
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
    if strcmp(method{nmethod},'three shuffles')
        load log_odd_wcorr
        log_odd_compare{7}{1} = log_odd;
        load log_odd_wcorr_global_remapped
        log_odd_compare{7}{2} = log_odd;

    elseif strcmp(method{nmethod},'two PRE shuffles')
        load log_odd_wcorr_two_PRE_shuffles
        log_odd_compare{6}{1} = log_odd;
        load log_odd_wcorr_two_PRE_shuffles_global_remapped
        log_odd_compare{6}{2} = log_odd;

    elseif strcmp(method{nmethod},'PRE place POST place')
        load log_odd_wcorr_PRE_place_POST_place
        log_odd_compare{5}{1} = log_odd;
        load log_odd_wcorr_PRE_place_POST_place_global_remapped
        log_odd_compare{5}{2} = log_odd;

    elseif strcmp(method{nmethod},'PRE place POST time')
        load log_odd_wcorr_PRE_place_POST_time
        log_odd_compare{4}{1} = log_odd;
        load log_odd_wcorr_PRE_place_POST_time_global_remapped
        log_odd_compare{4}{2} = log_odd;

    elseif strcmp(method{nmethod},'two POST')
        load log_odd_wcorr_two_POST_shuffles
        log_odd_compare{3}{1} = log_odd;
        load log_odd_wcorr_two_POST_shuffles_global_remapped
        log_odd_compare{3}{2} = log_odd;

    elseif strcmp(method{nmethod},'wcorr 1 shuffle + jump distance')
        load log_odd_wcorr_place_bin_circular_shift
        log_odd_compare{2}{1} = log_odd;
        load log_odd_wcorr_place_bin_circular_shift_global_remapped
        log_odd_compare{2}{2} = log_odd;

    elseif strcmp(method{nmethod},'place bin circular shift')
        load log_odd_wcorr_place_bin_circular_shift
        log_odd_compare{1}{1} = log_odd;
        load log_odd_wcorr_place_bin_circular_shift_global_remapped
        log_odd_compare{1}{2} = log_odd;

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
        data{nmethod}{2} = log_odd_compare{nmethod}{2}.normal_zscored.global_remapped_original

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
% alpha_level = linspace(0.1,0.75,length(low_threshold));

% log_odd_difference = [];
% log_odd_difference_CI = [];
% percent_sig_events = [];
% multi_event_percent = [];


cd ground_truth_original
if exist('log_odd_difference_multiple_shuffles.mat', 'file') ~= 2

    load log_odd_difference_single_shuffle
    % copy place bin shuffle data
    log_odd_difference1{1} = log_odd_difference{3};
    log_odd_difference_CI1{1} = log_odd_difference_CI{3};
    percent_sig_events1{1} = percent_sig_events{3};
    percent_sig_events_CI1{1} = percent_sig_events_CI{3};
    percent_multi_events1{1} = percent_multi_events{3};
    percent_multi_events_CI1{1} = percent_multi_events_CI{3};
    percent_shuffle_events1{1} = percent_shuffle_events{3};
    percent_shuffle_events_CI1{1} = percent_shuffle_events_CI{3};
    replay_score_compare1{1} = replay_score_compare{3};
    replay_score_compare_CI1{1} = replay_score_compare_CI{3};

    load log_odd_difference_jump_distance
    % copy wcorr 1 shuffle + jump distance 0.4 (From 71 because it equals to p <= 0.02)
    for epoch = 1:3
        for nshuffle = 1:2
            log_odd_difference1{2}{nshuffle}{epoch} = log_odd_difference{3}{nshuffle}{epoch}(71:end,:);
            log_odd_difference_CI1{2}{nshuffle}{epoch} = log_odd_difference_CI{3}{nshuffle}{epoch}(71:end,:);
            percent_sig_events1{2}{nshuffle}{epoch} = percent_sig_events{3}{nshuffle}{epoch}(71:end,:);
            percent_sig_events_CI1{2}{nshuffle}{epoch} = percent_sig_events_CI{3}{nshuffle}{epoch}(71:end,:);
            percent_multi_events1{2}{nshuffle}{epoch} = percent_multi_events{3}{nshuffle}{epoch}(71:end,:);
            percent_multi_events_CI1{2}{nshuffle}{epoch} = percent_multi_events_CI{3}{nshuffle}{epoch}(71:end,:);
            replay_score_compare1{2}{nshuffle}{epoch} = replay_score_compare{3}{nshuffle}{epoch}(71:end,:);
            replay_score_compare_CI1{2}{nshuffle}{epoch} = replay_score_compare_CI{3}{nshuffle}{epoch}(71:end,:);
        end
        percent_shuffle_events1{2}{epoch} = percent_shuffle_events{3}{epoch}(71:end,:);
        percent_shuffle_events_CI1{2}{epoch} = percent_shuffle_events_CI{3}{epoch}(71:end,:);
    end

    for nmethod = 3:length(method)
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

    save log_odd_difference_multiple_shuffles log_odd_difference log_odd_difference_CI percent_sig_events percent_sig_events_CI...
         percent_multi_events percent_multi_events_CI percent_shuffle_events percent_shuffle_events_CI replay_score_compare replay_score_compare_CI
    clear log_odd_difference1 log_odd_difference_CI1 percent_sig_events1 percent_sig_events_CI1 ...
        percent_multi_events1 percent_multi_events_CI1 percent_shuffle_events1 percent_shuffle_events_CI1 replay_score_compare1 replay_score_compare_CI1
    clear boot_log_odd_difference boot_percent_sig_events boot_percent_multi_events boot_percent_shuffle_events boot_replay_score
    cd ..
else
    load log_odd_difference_multiple_shuffles
    cd ..
end


cd ground_truth_original
if exist('log_odd_difference_multiple_shuffles_original.mat', 'file') ~= 2
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

                    if nmethod == 2 % Apply jump distance threshold 0.4
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

    save log_odd_difference_multiple_shuffles_original log_odd_difference_original percent_sig_events_original...
        percent_multi_events_original percent_shuffle_events_original replay_score_original
    clear boot_log_odd_difference boot_percent_sig_events boot_percent_multi_events boot_replay_score
    cd ..

else
    load log_odd_difference_multiple_shuffles_original
    cd ..
end

%% Demonstrate replay cross-validation framework

colour_line= {[8,81,156]/255,[67,147,195]/255,[244,109,67]/255,[215,48,39]/255};

fig = figure(1)
fig.Position = [834 116 850 700];

p_val_threshold = round(10.^[-3:0.01:log10(0.2)],4);
% p_val_threshold = 0.001:0.0001:0.2;
p_val_threshold = flip(unique(p_val_threshold));
p_val_threshold(1) = 0.2;
low_threshold = find(ismembertol(p_val_threshold,[0.0501 0.02 0.01 0.005 0.002 0.001],eps) == 1);
alpha_level = linspace(0.1,0.7,length(low_threshold));
size_level = linspace(40,80,length(low_threshold));

subplot(2,2,1)
for epoch = 2:3
    for nmethod = 4
        x = mean(log_odd_difference{nmethod}{1}{epoch},2)';
        y = mean(percent_sig_events{nmethod}{1}{epoch},2)';

        for threshold = 1:length(low_threshold)
            scatter(x(low_threshold(threshold)),y(low_threshold(threshold)),size_level(threshold),'filled','MarkerFaceColor',colour_line{epoch},'MarkerFaceAlpha',alpha_level(threshold),'MarkerEdgeColor',colour_line{epoch})
            hold on
        end

        y(isnan(x)) = [];
        x(isnan(x)) = [];
        UCI = log_odd_difference_CI{nmethod}{1}{epoch}(:,2)';
        UCI(isnan(UCI)) = [];
        LCI = log_odd_difference_CI{nmethod}{1}{epoch}(:,1)';
        LCI(isnan(LCI)) = [];
        p(epoch) = patch([LCI fliplr(UCI)], [y fliplr(y)], colour_line{epoch},'FaceAlpha','0.3','LineStyle','none');
        hold on

        xlabel('mean log odds difference')
        ylabel('Proportion of significant replay events')
    end

    ax = gca;
    ax.FontSize = 12;
    set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
    ylim([0 0.63])
    yticks([0:0.2:0.6 ])
end
legend([p(2),p(3)], {'RUN','POST'})


marker_shape = {'o','^','s'};
subplot(2,2,2)
for nmethod = 4
    for epoch = 2:3
        y = mean(percent_shuffle_events{nmethod}{epoch}(low_threshold,:),2);

        s(epoch) = scatter(p_val_threshold(low_threshold),y,20,marker_shape{epoch},...
            'MarkerFaceColor',colour_line{epoch},'MarkerEdgeColor',colour_line{epoch},'MarkerFaceAlpha','0.5')
        hold on
        m{epoch} = errorbar(p_val_threshold(low_threshold),y,y-percent_shuffle_events_CI{nmethod}{epoch}(low_threshold,1),...
            percent_shuffle_events_CI{nmethod}{epoch}(low_threshold,2)-y,'.','MarkerFaceColor',colour_line{epoch},'MarkerEdgeColor',colour_line{epoch})
        m{epoch}.Color = colour_line{epoch};
        %             hold on
        %             r(nshuffle) = rectangle('Position',pos)

        xlim([0 0.06])
        ylim([0 0.1])
        set(gca,'xscale','log')
        xticks(flip(p_val_threshold(low_threshold)))
        xticklabels({'0.001','0.002','0.005','0.01','0.02','0.05'})
        ylabel('Proportion of cell-id shuffled sig events')
        hold on
        xlabel('Sequenceness p value')
        %             title('Multitrack event proportion vs p value')
    end
    ax = gca;
    ax.FontSize = 12;
    set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
end

% Sig proportion and log odd at p value 0.05
subplot(2,2,3)
for nmethod = 4
    for epoch = 2:3
        x = mean(mean(log_odd_difference{nmethod}{1}{epoch}(low_threshold(1),:),2),2);
        y = mean(percent_sig_events{nmethod}{1}{epoch}(low_threshold(1),:),2);
        s(epoch) = scatter(x,y,20,marker_shape{epoch},...
            'MarkerFaceColor',colour_line{epoch},'MarkerEdgeColor',colour_line{epoch},'MarkerFaceAlpha','0.5')
        %             hold on
        %             m{nshuffle} = errorbar(p_val_threshold(low_threshold),y,y-percent_multi_events_CI{nmethod}{nshuffle}{epoch}(low_threshold,1),...
        %                 percent_multi_events_CI{nmethod}{nshuffle}{epoch}(low_threshold,2)-y,'.','MarkerFaceColor',colour_line{nshuffle},'MarkerEdgeColor',colour_line{nshuffle})
        %             m{nshuffle}.Color = colour_line{nshuffle};
        hold on
        error_x = [log_odd_difference_CI{nmethod}{1}{epoch}(low_threshold(1),1) log_odd_difference_CI{nmethod}{1}{epoch}(low_threshold(1),2)...
            log_odd_difference_CI{nmethod}{1}{epoch}(low_threshold(1),2) log_odd_difference_CI{nmethod}{1}{epoch}(low_threshold(1),1)]
        error_y = [percent_sig_events_CI{nmethod}{1}{epoch}(low_threshold(1),1) percent_sig_events_CI{nmethod}{1}{epoch}(low_threshold(1),1)...
            percent_sig_events_CI{nmethod}{1}{epoch}(low_threshold(1),2) percent_sig_events_CI{nmethod}{1}{epoch}(low_threshold(1),2)]

        patch(error_x,error_y,colour_line{epoch},'EdgeColor',colour_line{epoch},'FaceAlpha','0.2','EdgeAlpha','0.2')
        text{epoch} = sprintf('%s p =< %.3f (proportion = %.3f)',Behavioural_epoches{epoch},p_val_threshold(low_threshold(1)),mean(percent_shuffle_events{nmethod}{epoch}(low_threshold(1),:)));

        ylabel('Proportion of significant events')
        hold on
        xlabel('mean log odds difference')
        %             title(Behavioural_epoches{epoch})
        title('Original p value =< 0.05')
    end
end

% Sig proportion and log odd at shuffle-corrected p value 0.05
%         tempt = find(mean(percent_multi_events{nmethod}{nshuffle}{epoch},2) <= 0.05);
for nmethod = 4
    for epoch = 2:3
        [c index(epoch)] = min(abs(mean(percent_shuffle_events{nmethod}{epoch},2) - 0.05));
        x = mean(mean(log_odd_difference{nmethod}{1}{epoch}(index(epoch),:),2),2);
        y = mean(percent_sig_events{nmethod}{1}{epoch}(index(epoch),:),2);
        sc(epoch) = scatter(x,y,20,marker_shape{epoch},...
            'MarkerFaceColor',colour_line{epoch},'MarkerEdgeColor','k','MarkerFaceAlpha','0.5')
        %             hold on
        %             m{nshuffle} = errorbar(p_val_threshold(low_threshold),y,y-percent_multi_events_CI{nmethod}{nshuffle}{epoch}(low_threshold,1),...
        %                 percent_multi_events_CI{nmethod}{nshuffle}{epoch}(low_threshold,2)-y,'.','MarkerFaceColor',colour_line{nshuffle},'MarkerEdgeColor',colour_line{nshuffle})
        %             m{nshuffle}.Color = colour_line{nshuffle};
        hold on
        error_x = [log_odd_difference_CI{nmethod}{1}{epoch}(index(epoch),1) log_odd_difference_CI{nmethod}{1}{epoch}(index(epoch),2)...
            log_odd_difference_CI{nmethod}{1}{epoch}(index(epoch),2) log_odd_difference_CI{nmethod}{1}{epoch}(index(epoch),1)]
        error_y = [percent_sig_events_CI{nmethod}{1}{epoch}(index(epoch),1) percent_sig_events_CI{nmethod}{1}{epoch}(index(epoch),1)...
            percent_sig_events_CI{nmethod}{1}{epoch}(index(epoch),2) percent_sig_events_CI{nmethod}{1}{epoch}(index(epoch),2)]

        patch(error_x,error_y,colour_line{epoch},'EdgeColor',colour_line{epoch},'FaceAlpha','0.2','EdgeAlpha','0.7')
        textc{epoch} = sprintf('%s p =< %.3f (proportion = %.3f)',Behavioural_epoches{epoch},p_val_threshold(index(epoch)),mean(percent_shuffle_events{nmethod}{epoch}(index(epoch),:)));

        ylabel('Proportion of significant events')
        hold on
        xlabel('mean log odds difference')
        %             title(Behavioural_epoches{epoch})
        title('Original vs cell id shuffle corrected p value')
    end
end
%         plot([0 0.1],[0 0.3],'r')

xlim([0.5 2.2])
ylim([0 0.45])
ax = gca;
ax.FontSize = 12;
set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
yline(0.1,'--','LineWidth',0.01);
legend([s(2),s(3),sc(2),sc(3)], {text{2},text{3},textc{2},textc{3}})
cd ground_truth_original\Figure
filename = 'log odds framework original vs shuffle-corrected p value.pdf';
saveas(gcf,filename)
filename = 'log odds framework original vs shuffle-corrected p value.fig';
saveas(gcf,filename)
cd ..
cd ..
clf


%% 1 vs jump distance vs 2 vs 3 shuffles

shuffle_type = {'POST Place','POST Place + jump distance','POST place + POST time','PRE place + POST time','PRE place + POST place','PRE place + PRE spike','PRE place + PRE spike + POST Place'}';
p_val_threshold = round(10.^[-3:0.01:log10(0.2)],4);
% p_val_threshold = 0.001:0.0001:0.2;
p_val_threshold = flip(unique(p_val_threshold));
p_val_threshold(1) = 0.2;
low_threshold = find(ismembertol(p_val_threshold,[0.0501 0.02 0.01 0.005 0.002 0.001],eps) == 1);
colour_line= {[127,205,187]/255,[65,182,196]/255,'',[34,94,168]/255,'','',[37,52,148]/255};
alpha_level = linspace(0.2,0.7,length(low_threshold));
size_level = linspace(40,80,length(low_threshold));

fig = figure(2)
fig.Position = [834 116 850 700];
nfig = 1;
for epoch = 1:3
    subplot(2,2,nfig)
    for nmethod = [1 2 4 7]
        x = mean(log_odd_difference{nmethod}{1}{epoch},2)';
        y = mean(percent_sig_events{nmethod}{1}{epoch},2)';

        for threshold = 1:length(low_threshold)
            scatter(x(low_threshold(threshold)),y(low_threshold(threshold)),size_level(threshold),'filled','MarkerFaceColor',colour_line{nmethod},'MarkerFaceAlpha',alpha_level(threshold),'MarkerEdgeColor',colour_line{nmethod})
            hold on
        end

        UCI = log_odd_difference_CI{nmethod}{1}{epoch}(:,2)';
        UCI(isnan(x)) = [];
        LCI = log_odd_difference_CI{nmethod}{1}{epoch}(:,1)';
        LCI(isnan(x)) = [];
        y(isnan(x)) = [];
        x(isnan(x)) = [];
        p(nmethod) = patch([LCI fliplr(UCI)], [y fliplr(y)], colour_line{nmethod},'FaceAlpha','0.3','LineStyle','none');
        hold on

        xlabel('mean log odds difference')
        ylabel('Proportion of significant replay events')
  
    title(Behavioural_epoches{epoch});
    end
    ylim([0 0.7])
    yticks([0 0.2 0.4 0.6])
    ax = gca;
    ax.FontSize = 12;
    set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
    nfig = nfig + 1;

end
%     legend([p(1) p(3) p(6)], {'Shuffle 3','Shuffle 2+3','Shuffle 1+2+3'})
legend([p(1),p(2),p(4),p(7)], {shuffle_type{1},shuffle_type{2},shuffle_type{4},shuffle_type{7}},'Position',[0.7 0.3 0.05 0.05])


cd ground_truth_original\Figure
filename = sprintf('wcorr 3 shuffle comparisions CI.pdf')
saveas(gcf,filename)
filename = sprintf('wcorr 3 shuffle comparisions CI.fig')
saveas(gcf,filename)
cd ..
cd ..
clf
%


p_val_threshold = round(10.^[-3:0.01:log10(0.2)],4);
p_val_threshold = flip(unique(p_val_threshold));
p_val_threshold(1) = 0.2;
low_threshold = find(ismembertol(p_val_threshold,[0.0501 0.02 0.01 0.005 0.002 0.001],eps) == 1);
colour_line= {[127,205,187]/255,[65,182,196]/255,'',[34,94,168]/255,'','',[37,52,148]/255};
Behavioural_epoches = {'PRE','RUN','POST'};
marker_shape = {'o','^','s'};

fig = figure(2)
fig.Position = [834 116 650 650];
nfig = 1;
for epoch = 2:3
    for nmethod = [1 2 4 7]

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
    ax = gca;
    ax.FontSize = 12;
    set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
    nfig = nfig + 1;
end
legend([s(1),s(2),s(4),s(7)], {shuffle_type{1},shuffle_type{2},shuffle_type{4},shuffle_type{7}},'Position',[0.7 0.3 0.05 0.05])

%     legend([s(1),s(2),s(3),s(4)], {'spike train circular shift','place field circular shift','place bin circular shift','time bin permutation'},'Position',[0.35 0.2 0.05 0.05])
%     sgtitle('Multitrack event proportion vs Sequenceness p value')
cd ground_truth_original\Figure
filename = sprintf('wcorr 3 shuffles cell-id shuffled sig proportion.pdf')
saveas(gcf,filename)
filename = sprintf('wcorr 3 shuffles cell-id shuffled sig proportion.fig')
saveas(gcf,filename)
cd ..
cd ..
clf


% Inset for false positive rate at p value 0.05 and p value for false
% positive rate 0.05
fig = figure(3)
fig.Position = [834 116 280 470];
count = 1;

subplot(2,2,1)
for nmethod = [1 2 4 7]
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
ylim([0 0.25])
xlim([0 2])
set(gca,'xtick',[]);

subplot(2,2,2)
count = 1;
for nmethod = [1 2 4 7]
    for epoch = 2:3
        [xx index] = min(abs(mean(percent_shuffle_events{nmethod}{epoch},2) - 0.05));
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
ylim([0 0.2])
xlim([0 2])
set(gca,'xtick',[]);
cd ground_truth_original\Figure
filename = sprintf('wcorr 3 shuffles false positive rate and p value inset (POST ripple).pdf')
saveas(gcf,filename)
filename = sprintf('wcorr 3 shuffles false positive rate and p value inset (POST ripple).fig')
saveas(gcf,filename)
cd ..
cd ..
clf


% Sig proportion and log odd at p value 0.05
fig = figure(3)
fig.Position = [834 116 850 885];
subplot(2,2,1)
count = 1;
for nmethod = [1 2 4 7]
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
        text{count} = sprintf('%s p =< %.3f (proportion = %.3f)',shuffle_type{nmethod},p_val_threshold(low_threshold(1)),mean(percent_shuffle_events{nmethod}{epoch}(low_threshold(1),:)));

        ylabel('Proportion of significant events')
        hold on
        xlabel('mean log odds difference')
        %             title(Behavioural_epoches{epoch})
        title('Original p value =< 0.05')
        count = count + 1;
    end
end
xlim([0.5 3])
ylim([0 0.55])
yticks([0.1 0.2 0.3 0.4 0.5])
ax = gca;
ax.FontSize = 12;
set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
yline(0.1,'--')
legend([s(1:8)], {text{1:8}})


% Sig proportion and log odd at shuffle-corrected p value 0.05
subplot(2,2,2)
count = 1;
%         tempt = find(mean(percent_multi_events{nmethod}{nshuffle}{epoch},2) <= 0.05);
for nmethod = [1 2 4 7]
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
        text{count} = sprintf('%s p =< %.3f (proportion = %.3f)',shuffle_type{nmethod},p_val_threshold(index(nmethod)),mean(percent_shuffle_events{nmethod}{epoch}(index(nmethod),:)));
%         text{nmethod} = sprintf('%s p =< %.3f (proportion = %.3f)',shuffle_type{nmethod},p_val_threshold(index(nmethod)),mean(percent_multi_events{nmethod}{1}{epoch}(index(nmethod),:)));

        ylabel('Proportion of significant events')
        hold on
        xlabel('mean log odds difference')
        %             title(Behavioural_epoches{epoch})
        title('Equivalent p value when cell-id shuffled sig proportion = 0.05')
        count = count + 1;
    end
end
%         plot([0 0.1],[0 0.3],'r')
xlim([0.5 3])
ylim([0 0.4])
yticks([0.1 0.2 0.3 0.4])
ax = gca;
ax.FontSize = 12;
set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
yline(0.1,'--')
legend([s(1:8)], {text{1:8}})

cd ground_truth_original\Figure
filename = sprintf('wcorr 3 shuffles original vs shuffle-corrected p value.pdf')
saveas(gcf,filename)
cd ..
cd ..
clf


% 
% %% 2 shuffle combinations
% 
% fig = figure(2)
% fig.Position = [834 116 850 700];
% 
% shuffle_type = {'POST Place','POST Place + jump distance','POST place + POST time','PRE place + POST time','PRE place + POST place','PRE place + PRE spike','PRE place + PRE spike + POST Place'}';
% p_val_threshold = round(10.^[-3:0.01:log10(0.2)],4);
% p_val_threshold = flip(unique(p_val_threshold));
% p_val_threshold(1) = 0.2;
% low_threshold = find(ismembertol(p_val_threshold,[0.0501 0.02 0.01 0.005 0.002 0.001],eps) == 1);
% colour_line = {'','',[127,205,187]/255,[34,94,168]/255,[65,182,196]/255,[37,52,148]/255,''};
% alpha_level = linspace(0.1,0.7,length(low_threshold));
% %     colour_line = {[37,52,148]/255,[34,94,168]/255,[29,145,192]/255,[65,182,196]/255};
% %     colour_line = flip({[12,44,132]/255,'',[29,145,192]/255,'','',[127,205,187]/255});
% nfig = 1;
% for epoch = 2:3
%     subplot(2,2,nfig)
%     for nmethod = 3:6
%         x = mean(log_odd_difference{nmethod}{1}{epoch},2)';
%         y = mean(percent_sig_events{nmethod}{1}{epoch},2)';
% 
%         for threshold = 1:length(low_threshold)
%             scatter(x(low_threshold(threshold)),y(low_threshold(threshold)),'filled','MarkerFaceColor',colour_line{nmethod},...
%                 'MarkerFaceAlpha',alpha_level(threshold),'MarkerEdgeColor',colour_line{nmethod})
%             hold on
%         end
% 
%         UCI = log_odd_difference_CI{nmethod}{1}{epoch}(:,2)';
%         UCI(isnan(x)) = [];
%         LCI = log_odd_difference_CI{nmethod}{1}{epoch}(:,1)';
%         LCI(isnan(x)) = [];
% 
%         y(isnan(x)) = [];
%         x(isnan(x)) = [];
%         p(nmethod) = patch([LCI fliplr(UCI)], [y fliplr(y)], colour_line{nmethod},'FaceAlpha','0.3','LineStyle','none');
%         hold on
%         %             plot(x,y,colour_line{nshuffle});
% 
%         xlabel('mean log odds difference')
%         ylabel('Proportion of significant replay events')
%     end
% 
%     ylim([0 0.7])
%     title(Behavioural_epoches{epoch});
%     nfig = nfig + 1;
% end
% 
% %     legend([p(1) p(3) p(6)], {'Shuffle 3','Shuffle 2+3','Shuffle 1+2+3'})
% legend([p(3) p(5) p(4) p(6)], {shuffle_type{3},shuffle_type{5},shuffle_type{4},shuffle_type{6}},'Position',[0.7 0.3 0.05 0.05])
% 
% %     legend([p(1) p(3) p(6)], {'PRE Place','PRE Place + POST Place','PRE Place + PRE Spike + POST Place'})
% %     legend([p(1),p(2),p(3),p(4),p(5),p(6)], {'Shuffle 3','Shuffle 3+4','Shuffle 2+3','Shuffle 1+3','Shuffle 1+2','Shuffle 1+2+3'})
% cd ground_truth_original\Figure
% filename = '2 shuffles log odds comparision CI.pdf'
% saveas(gcf,filename)
% filename = '2 shuffles log odds comparision CI.fig'
% saveas(gcf,filename)
% cd ..
% cd ..
% clf
% %
% 
% 
% 
% %     colour_line = {'',[127,205,187]/255,[34,94,168]/255,[29,145,192]/255,[37,52,148]/255,''};
% % colour_line = {'',[127,205,187]/255,[65,182,196]/255,[34,94,168]/255,[37,52,148]/255,''};
% colour_line = {'','',[127,205,187]/255,[34,94,168]/255,[65,182,196]/255,[37,52,148]/255,''};
% p_val_threshold = round(10.^[-3:0.01:log10(0.2)],4);
% p_val_threshold = flip(unique(p_val_threshold));
% p_val_threshold(1) = 0.2;
% low_threshold = find(ismembertol(p_val_threshold,[0.0501 0.02 0.01 0.005 0.002 0.001],eps) == 1);
% shuffle_type = {'POST Place','wcorr 1 shuffle + jump distance','POST place + POST time','PRE place + POST time','PRE place + POST place','PRE place + PRE spike','PRE place + PRE spike + POST Place'}';
% Behavioural_epoches = {'PRE','RUN','POST'};
% marker_shape = {'o','^','s'};
% nfig = 1;
% fig = figure(2)
% fig.Position = [834 116 850 885];
% 
% for epoch = 2:3
%     for nmethod = 3:6
%         subplot(2,2,nfig)
%         y = mean(percent_shuffle_events{nmethod}{epoch}(low_threshold,:),2);
% 
%         s(nmethod) = scatter(p_val_threshold(low_threshold),y,20,marker_shape{epoch},...
%             'MarkerFaceColor',colour_line{nmethod},'MarkerEdgeColor',colour_line{nmethod},'MarkerFaceAlpha','0.5')
%         hold on
%         m{nmethod} = errorbar(p_val_threshold(low_threshold),y,y-percent_shuffle_events_CI{nmethod}{epoch}(low_threshold,1),...
%             percent_shuffle_events_CI{nmethod}{epoch}(low_threshold,2)-y,'.','MarkerFaceColor',colour_line{nmethod},'MarkerEdgeColor',colour_line{nmethod})
%         m{nmethod}.Color = colour_line{nmethod};
%         %             hold on
%         %             r(nshuffle) = rectangle('Position',pos)
%         set(gca,'xscale','log')
%         xticks(flip(p_val_threshold(low_threshold)))
%         xticklabels({'0.001','0.002','0.005','0.01','0.02','0.05'})
%         xlim([0 0.06])
%         ylim([0 0.2])
%         ylabel('Proportion of cell-id shuffled events')
%         hold on
%         xlabel('Sequenceness p value')
%         title(Behavioural_epoches{epoch})
%         %             title('Multitrack event proportion vs p value')
%     end
%     nfig = nfig + 1;
% end
% legend([s(3) s(5) s(4) s(6)], {shuffle_type{3},shuffle_type{5},shuffle_type{4},shuffle_type{6}})
% 
% %     legend([s(1),s(2),s(3),s(4)], {'spike train circular shift','place field circular shift','place bin circular shift','time bin permutation'},'Position',[0.35 0.2 0.05 0.05])
% %     sgtitle('Multitrack event proportion vs Sequenceness p value')
% cd ground_truth_original\Figure
% filename = sprintf('wcorr 2 shuffles cell-id shuffled sig proportion.pdf')
% saveas(gcf,filename)
% filename = sprintf('wcorr 2 shuffles cell-id shuffled sig proportion.fig')
% saveas(gcf,filename)
% cd ..
% cd ..
% clf
% 
% 
% 
% % Sig proportion and log odd at p value 0.05
% % colour_line = {'',[127,205,187]/255,[65,182,196]/255,[34,94,168]/255,[37,52,148]/255,''};
% colour_line = {'','',[127,205,187]/255,[34,94,168]/255,[65,182,196]/255,[37,52,148]/255,''};
% shuffle_type = {'POST Place','wcorr 1 shuffle + jump distance','POST place + POST time','PRE place + POST time','PRE place + POST place','PRE place + PRE spike','PRE place + PRE spike + POST Place'}';
% Behavioural_epoches = {'PRE','RUN','POST'};
% marker_shape = {'o','^','s'};
% 
% fig = figure(3)
% fig.Position = [834 116 850 885];
% subplot(2,2,1)
% count = 1;
% for nmethod = 3:6
%     for epoch = 2:3
%         x = mean(mean(log_odd_difference{nmethod}{1}{epoch}(low_threshold(1),:),2),2);
%         y = mean(percent_sig_events{nmethod}{1}{epoch}(low_threshold(1),:),2);
%         s(count) = scatter(x,y,20,marker_shape{epoch},...
%             'MarkerFaceColor',colour_line{nmethod},'MarkerEdgeColor',colour_line{nmethod},'MarkerFaceAlpha','0.5')
%         %             hold on
%         %             m{nshuffle} = errorbar(p_val_threshold(low_threshold),y,y-percent_multi_events_CI{nmethod}{nshuffle}{epoch}(low_threshold,1),...
%         %                 percent_multi_events_CI{nmethod}{nshuffle}{epoch}(low_threshold,2)-y,'.','MarkerFaceColor',colour_line{nshuffle},'MarkerEdgeColor',colour_line{nshuffle})
%         %             m{nshuffle}.Color = colour_line{nshuffle};
%         hold on
%         error_x = [log_odd_difference_CI{nmethod}{1}{epoch}(low_threshold(1),1) log_odd_difference_CI{nmethod}{1}{epoch}(low_threshold(1),2)...
%             log_odd_difference_CI{nmethod}{1}{epoch}(low_threshold(1),2) log_odd_difference_CI{nmethod}{1}{epoch}(low_threshold(1),1)]
%         error_y = [percent_sig_events_CI{nmethod}{1}{epoch}(low_threshold(1),1) percent_sig_events_CI{nmethod}{1}{epoch}(low_threshold(1),1)...
%             percent_sig_events_CI{nmethod}{1}{epoch}(low_threshold(1),2) percent_sig_events_CI{nmethod}{1}{epoch}(low_threshold(1),2)]
% 
%         patch(error_x,error_y,colour_line{nmethod},'EdgeColor',colour_line{nmethod},'FaceAlpha','0.2','EdgeAlpha','0.2')
%         text{count} = sprintf('%s p =< %.3f (proportion = %.3f)',shuffle_type{nmethod},p_val_threshold(low_threshold(1)),mean(percent_shuffle_events{nmethod}{epoch}(low_threshold(1),:)));
% 
%         ylabel('Proportion of significant events')
%         hold on
%         xlabel('mean log odds difference')
%         %             title(Behavioural_epoches{epoch})
%         title('Original p value =< 0.05')
% 
%         count = count +1;
%     end
% end
% legend([s(1:8)], {text{1:8}})
% 
% xlim([-1 3])
% ylim([0 0.7])
% 
% % Sig proportion and log odd at shuffle-corrected p value 0.05
% subplot(2,2,2)
% count = 1;
% %         tempt = find(mean(percent_multi_events{nmethod}{nshuffle}{epoch},2) <= 0.05);
% for nmethod = 3:6
%     for epoch = 2:3
%         [c index(nmethod)] = min(abs(mean(percent_shuffle_events{nmethod}{epoch},2) - 0.05));
%         x = mean(mean(log_odd_difference{nmethod}{1}{epoch}(index(nmethod),:),2),2);
%         y = mean(percent_sig_events{nmethod}{1}{epoch}(index(nmethod),:),2);
%         s(count) = scatter(x,y,20,marker_shape{epoch},...
%             'MarkerFaceColor',colour_line{nmethod},'MarkerEdgeColor',colour_line{nmethod},'MarkerFaceAlpha','0.5')
%         %             hold on
%         %             m{nshuffle} = errorbar(p_val_threshold(low_threshold),y,y-percent_multi_events_CI{nmethod}{nshuffle}{epoch}(low_threshold,1),...
%         %                 percent_multi_events_CI{nmethod}{nshuffle}{epoch}(low_threshold,2)-y,'.','MarkerFaceColor',colour_line{nshuffle},'MarkerEdgeColor',colour_line{nshuffle})
%         %             m{nshuffle}.Color = colour_line{nshuffle};
%         hold on
%         error_x = [log_odd_difference_CI{nmethod}{1}{epoch}(index(nmethod),1) log_odd_difference_CI{nmethod}{1}{epoch}(index(nmethod),2)...
%             log_odd_difference_CI{nmethod}{1}{epoch}(index(nmethod),2) log_odd_difference_CI{nmethod}{1}{epoch}(index(nmethod),1)]
%         error_y = [percent_sig_events_CI{nmethod}{1}{epoch}(index(nmethod),1) percent_sig_events_CI{nmethod}{1}{epoch}(index(nmethod),1)...
%             percent_sig_events_CI{nmethod}{1}{epoch}(index(nmethod),2) percent_sig_events_CI{nmethod}{1}{epoch}(index(nmethod),2)]
% 
%         patch(error_x,error_y,colour_line{nmethod},'EdgeColor',colour_line{nmethod},'FaceAlpha','0.2','EdgeAlpha','0.2')
%         text{nmethod} = sprintf('%s p =< %.3f (proportion = %.3f)',shuffle_type{nmethod},p_val_threshold(index(nmethod)),mean(percent_shuffle_events{nmethod}{epoch}(index(nmethod),:)));
% 
%         ylabel('Proportion of significant events')
%         hold on
%         xlabel('mean log odds difference')
%         %             title(Behavioural_epoches{epoch})
%         title('Equivalent p value when cell-id shuffled sig proportion = 0.05')
%         count = count + 1;
%     end
% end
% %         plot([0 0.1],[0 0.3],'r')
% legend([s(1:8)], {text{1:8}})
% 
% xlim([-1 3])
% ylim([0 0.7])
% 
% cd ground_truth_original\Figure
% filename = sprintf('wcorr 2 shuffles original vs shuffle-corrected p value.pdf')
% saveas(gcf,filename)
% cd ..
% cd ..
% clf

%% Wcorr 1 vs 2 vs 3 shuffles table
%     shuffle_type = {'place field circular shift','spike train circular shift','place bin circular shift','time bin circular shift'}';
shuffle_type = {'place bin circular shift','wcorr 1 shuffle + jump distance','two POST','PRE place POST time','PRE place POST place','two PRE shuffles','three shuffles'};

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
shuffle_to_compare = [1 2 4 7];
for epoch = 1:3
    %     tempt = [mean(percent_sig_events{nmethod}{1}{epoch}(low_threshold(1),:)) mean(percent_sig_events{nmethod}{2}{epoch}(low_threshold(1),:))...
    %         mean(percent_sig_events{nmethod}{3}{epoch}(low_threshold(1),:)) mean(percent_sig_events{nmethod}{4}{epoch}(low_threshold(1),:))];
    %     [minimum mindex] = min(tempt);
    for nmethod = shuffle_to_compare
        for nshuffle = 1:1
            shuffle_name(c) = {shuffle_type{nmethod}};
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

writetable(Table,'1 vs 2 vs 3 shuffle combination effect table.xlsx')
cd ..
cd ..


%% Wcorr shuffle combination
shuffle_type = {'place bin circular shift','wcorr 1 shuffle + jump distance','two POST','PRE place POST time','PRE place POST place','two PRE shuffles','three shuffles'};

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
shuffle_to_compare = [3 5 4 6];
for epoch = 1:3
    %     tempt = [mean(percent_sig_events{nmethod}{1}{epoch}(low_threshold(1),:)) mean(percent_sig_events{nmethod}{2}{epoch}(low_threshold(1),:))...
    %         mean(percent_sig_events{nmethod}{3}{epoch}(low_threshold(1),:)) mean(percent_sig_events{nmethod}{4}{epoch}(low_threshold(1),:))];
    %     [minimum mindex] = min(tempt);
    for nmethod = shuffle_to_compare
        for nshuffle = 1:1
            shuffle_name(c) = {shuffle_type{nmethod}};
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

writetable(Table,'2 shuffle combination effect table.xlsx')
cd ..
cd ..



end


