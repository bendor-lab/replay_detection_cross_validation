function plot_ground_truth_cell_id_shuffle_validation(folders,option,method)
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

total_number{3} = total_number{2};
total_number{4} = total_number{2};
total_number{5} = total_number{2};

index = [];
epoch_index = [];

for nmethod = 1:length(method)

    cd ground_truth_original
    if strcmp(method{nmethod},'time shuffle')
        load log_odd_wcorr_time_bin_permutation
        log_odd_compare{4}{1} = log_odd;
        load log_odd_wcorr_time_bin_permutation_global_remapped
        log_odd_compare{4}{2} = log_odd;
        load log_odd_wcorr_time_bin_permutation_cross_experiment_shuffled
        log_odd_compare{4}{3} = log_odd;
        load log_odd_wcorr_time_bin_permutation_place_field_shifted
        log_odd_compare{4}{4} = log_odd;
        load log_odd_wcorr_time_bin_permutation_spike_train_shifted
        log_odd_compare{4}{5} = log_odd;

    elseif strcmp(method{nmethod},'place shuffle')
        load log_odd_wcorr_place_bin_circular_shift
        log_odd_compare{3}{1} = log_odd;
        load log_odd_wcorr_place_bin_circular_shift_global_remapped
        log_odd_compare{3}{2} = log_odd;
        load log_odd_wcorr_place_bin_circular_shift_cross_experiment_shuffled
        log_odd_compare{3}{3} = log_odd;
        load log_odd_wcorr_place_bin_circular_shift_place_field_shifted
        log_odd_compare{3}{4} = log_odd;
        load log_odd_wcorr_place_bin_circular_shift_spike_train_shifted
        log_odd_compare{3}{5} = log_odd;

    elseif strcmp(method{nmethod},'place field shuffle')
        load log_odd_wcorr_place_field_circular_shift
        log_odd_compare{2}{1} = log_odd;
        load log_odd_wcorr_place_field_circular_shift_global_remapped
        log_odd_compare{2}{2} = log_odd;
        load log_odd_wcorr_place_field_circular_shift_cross_experiment_shuffled
        log_odd_compare{2}{3} = log_odd;
        load log_odd_wcorr_place_field_circular_shift_place_field_shifted
        log_odd_compare{2}{4} = log_odd;
        load log_odd_wcorr_place_field_circular_shift_spike_train_shifted
        log_odd_compare{2}{5} = log_odd;

    elseif strcmp(method{nmethod},'spike shuffle')
        load log_odd_wcorr_spike_train_circular_shift
        log_odd_compare{1}{1} = log_odd;
        load log_odd_wcorr_spike_train_circular_shift_global_remapped
        log_odd_compare{1}{2} = log_odd;
        load log_odd_wcorr_spike_train_circular_shift_cross_experiment_shuffled
        log_odd_compare{1}{3} = log_odd;
        load log_odd_wcorr_spike_train_circular_shift_place_field_shifted
        log_odd_compare{1}{4} = log_odd;
        load log_odd_wcorr_spike_train_circular_shift_spike_train_shifted
        log_odd_compare{1}{5} = log_odd;


    elseif strcmp(method{nmethod},'cell id shuffle')
        load log_odd_wcorr_cell_id_shuffle
        log_odd_compare{5}{1} = log_odd;
        load log_odd_wcorr_cell_id_shuffle_global_remapped
        log_odd_compare{5}{2} = log_odd;
        load log_odd_wcorr_cell_id_shuffle_cross_experiment_shuffled
        log_odd_compare{5}{3} = log_odd;
        load log_odd_wcorr_cell_id_shuffle_place_field_shifted
        log_odd_compare{5}{4} = log_odd;
        load log_odd_wcorr_cell_id_shuffle_spike_train_shifted
        log_odd_compare{5}{5} = log_odd;

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
            replay_score{nmethod}{nshuffle} = log_odd_compare{nmethod}{nshuffle}.best_score;
        end

        %     elseif strcmp(option,'original')
        %
        %         data{nmethod}{1} = log_odd_compare{nmethod}{1}.normal_zscored.original;
        %         data{nmethod}{2} = log_odd_compare{nmethod}{2}.normal_zscored.global_remapped_original
        %
        %         for nshuffle = 1:length(log_odd_compare{nmethod})
        %             log_pval{nmethod}{nshuffle} = log_odd_compare{nmethod}{nshuffle}.pvalue;
        %             segment_id{nmethod}{nshuffle} = log_odd_compare{nmethod}{nshuffle}.segment_id;
        %             replay_score{nmethod}{nshuffle} = log_odd_compare{nmethod}{nshuffle}.best_score;
        %         end
    end

    %     colour_line2 = {'k--','b--','r--','g--'};
    colour_symbol={'bo','ro','go','ko'};
    Behavioural_epoches = {'PRE','RUN','POST'};

end

p_val_threshold = round(10.^[-3:0.01:log10(0.2)],4);
% p_val_threshold = 0.001:0.0001:0.2;
p_val_threshold = flip(unique(p_val_threshold));
p_val_threshold(1) = 0.2;

log_odd_difference = [];
log_odd_difference_CI = [];
percent_sig_events = [];

cd ground_truth_original
if exist('log_odd_difference_single_shuffle_validation.mat', 'file') ~= 2
    for nmethod = 1:length(method)
        for epoch = 1:3
            for nshuffle = 1:length(log_pval{nmethod})
                %                 multi_event_index = find(log_odd_compare{nshuffle}.multi_track == 1);
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

                        if nshuffle ~= 1
                            % cell id shuffled mean significant event proportion
                            boot_percent_shuffle_events(threshold,nboot) = ((length(track_1_index)+length(track_2_index))/2) / total_number{nshuffle}(epoch);
                        end
                    end
                end

                log_odd_difference{nmethod}{nshuffle}{epoch} = boot_log_odd_difference;
                log_odd_difference_CI{nmethod}{nshuffle}{epoch} = prctile(boot_log_odd_difference,[2.5 97.5],2);
                percent_sig_events{nmethod}{nshuffle}{epoch} = boot_percent_sig_events;
                percent_sig_events_CI{nmethod}{nshuffle}{epoch} = prctile(boot_percent_sig_events,[2.5 97.5],2);
                percent_multi_events{nmethod}{nshuffle}{epoch} = boot_percent_multi_events;
                percent_multi_events_CI{nmethod}{nshuffle}{epoch} = prctile(boot_percent_multi_events,[2.5 97.5],2);
                replay_score_compare{nmethod}{nshuffle}{epoch} = boot_replay_score;
                replay_score_compare_CI{nmethod}{nshuffle}{epoch} = prctile(boot_replay_score,[2.5 97.5],2);

                if nshuffle ~= 1
                    percent_shuffle_events{nmethod}{epoch}{nshuffle} = boot_percent_shuffle_events;
                    percent_shuffle_events_CI{nmethod}{epoch}{nshuffle} = prctile(boot_percent_shuffle_events,[2.5 97.5],2);
                end

                toc
            end
        end
    end

    save log_odd_difference_single_shuffle_validation log_odd_difference log_odd_difference_CI percent_sig_events percent_sig_events_CI...
        percent_multi_events percent_multi_events_CI percent_shuffle_events percent_shuffle_events_CI replay_score_compare replay_score_compare_CI
    clear boot_log_odd_difference boot_percent_sig_events boot_percent_multi_events boot_percent_shuffle_events boot_replay_score
    cd ..
else
    load log_odd_difference_single_shuffle_validation
    cd ..
end


log_odd_difference_original = [];
percent_sig_events_original = [];

cd ground_truth_original
if exist('log_odd_difference_single_shuffle_validation_original.mat', 'file') ~= 2
    for nmethod = 1:length(method)
        for epoch = 1:3
            for nshuffle = 1:length(log_pval{nmethod})
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

                    if nshuffle ~= 1
                        % cell id shuffled mean significant event proportion
                        boot_percent_shuffle_events(threshold) = ((length(track_1_index)+length(track_2_index))/2) / total_number{nshuffle}(epoch);
                    end
                end

                log_odd_difference_original{nmethod}{nshuffle}{epoch} = boot_log_odd_difference;
                percent_sig_events_original{nmethod}{nshuffle}{epoch} = boot_percent_sig_events;
                percent_multi_events_original{nmethod}{nshuffle}{epoch} = boot_percent_multi_events;
                replay_score_original{nmethod}{nshuffle}{epoch} = boot_replay_score;

                if nshuffle ~= 1
                    percent_shuffle_events_original{nmethod}{epoch} = boot_percent_shuffle_events;
                end

                toc
            end

        end
    end

    save log_odd_difference_single_shuffle_validation_original log_odd_difference_original percent_sig_events_original...
        percent_multi_events_original percent_shuffle_events_original replay_score_original
    clear boot_replay_score boot_percent_multi_events boot_percent_sig_events boot_log_odd_difference
    cd ..
else
    load log_odd_difference_single_shuffle_validation_original
    cd ..
end

p_val_threshold = round(10.^[-3:0.01:log10(0.2)],4);
% p_val_threshold = 0.001:0.0001:0.2;
p_val_threshold = flip(unique(p_val_threshold));
p_val_threshold(1) = 0.2;
low_threshold = find(ismembertol(p_val_threshold,[0.0501 0.02 0.01 0.005 0.002 0.001],eps) == 1);

colour_line= {[8,81,156]/255,[67,147,195]/255,[244,109,67]/255,[215,48,39]/255,[94,60,153]/255};
% colour_line= {[50,136,189]/255,[26,152,80]/255,[244,109,67]/255,[213,62,79]/255};
alpha_level = linspace(0.2,0.7,length(low_threshold));
size_level = linspace(40,80,length(low_threshold));
% size_level = linspace(40,40,length(low_threshold));

ndataset = [] %

shuffle_type = {'spike train shuffle','place field shuffle','place bin shuffle','time bin shuffle','cell id shuffle'}';
%
% for nmethod = [5 4 3 2 1]
%     fig = figure(nmethod)
%     fig.Position = [834 116 850 700];
%     nfig = 1;
%     for epoch = 1:3
%         subplot(2,2,nfig)
%         for ndataset = [1 2 3 4]
%             x = mean(log_odd_difference{nmethod}{ndataset}{epoch},2)';
%             y = mean(percent_sig_events{nmethod}{ndataset}{epoch},2)';
%
%             for threshold = 1:length(low_threshold)
%                 scatter(x(low_threshold(threshold)),y(low_threshold(threshold)),size_level(threshold),'filled','MarkerFaceColor',colour_line{ndataset},'MarkerFaceAlpha',alpha_level(threshold),'MarkerEdgeColor',colour_line{ndataset})
%                 hold on
%             end
%
%
%             UCI = log_odd_difference_CI{nmethod}{ndataset}{epoch}(:,2)';
%             UCI(isnan(x)) = [];
%             LCI = log_odd_difference_CI{nmethod}{ndataset}{epoch}(:,1)';
%             LCI(isnan(x)) = [];
%             y(isnan(x)) = [];
%             x(isnan(x)) = [];
%             p(ndataset) = patch([LCI fliplr(UCI)], [y fliplr(y)], colour_line{ndataset},'FaceAlpha','0.3','LineStyle','none');
%             hold on
%
%
%
%             xlabel('mean log odds difference')
%             ylabel('Proportion of significant replay events')
%         end
%
%         ylim([0 0.8])
%         title(Behavioural_epoches{epoch});
%         ax = gca;
%         ax.FontSize = 12;
%         set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
%         nfig = nfig + 1;
%     end
%     legend([p(1),p(2),p(3),p(4)], {'Original','Cell id randomised data',...
%         'cross experiment shuffled data','place field shuffled data'})
%     sgtitle(shuffle_type{nmethod})
% end
%
% % legend([p(1),p(2),p(3),p(4)], {'spike train circular shift','place field circular shift','place bin circular shift','time bin permutation'})
%
% cd ground_truth_original\Figure
% filename = 'wcorr shuffle validation log odds comparision CI.pdf';
% saveas(gcf,filename)
% filename = 'wcorr shuffle validation log odds comparision CI.fig';
% saveas(gcf,filename)
% cd ..
% cd ..
% clf


p_val_threshold = round(10.^[-3:0.01:log10(0.2)],4);
p_val_threshold = flip(unique(p_val_threshold));
p_val_threshold(1) = 0.2;
low_threshold = find(ismembertol(p_val_threshold,[0.0501 0.02 0.01 0.005 0.002 0.001],eps) == 1);

shuffle_type = {'spike train shuffle','place field shuffle','place bin shuffle','time bin shuffle','cell id shuffle'}';
dataset_type = {'original','cell id shuffled dataset','cross experiment shuffled dataset','place field shifted dataset','spike train shifted dataset'};

Behavioural_epoches = {'PRE','RUN','POST'};
marker_shape = {'o','^','s'};


for nmethod = 1:5
    fig = figure(nmethod)
    fig.Position = [834 116 850 885];
    nfig = 1;

    for epoch = 1:3
        for ndataset = [2 4 5]

            subplot(2,2,nfig)
            y = mean(percent_shuffle_events{nmethod}{epoch}{ndataset}(low_threshold,:),2);

            s(ndataset) = scatter(p_val_threshold(low_threshold),y,20,marker_shape{epoch},...
                'MarkerFaceColor',colour_line{ndataset},'MarkerEdgeColor',colour_line{ndataset},'MarkerFaceAlpha','0.5')
            hold on
            m{ndataset} = errorbar(p_val_threshold(low_threshold),y,y-percent_shuffle_events_CI{nmethod}{epoch}{ndataset}(low_threshold,1),...
                percent_shuffle_events_CI{nmethod}{epoch}{ndataset}(low_threshold,2)-y,'.','MarkerFaceColor',colour_line{ndataset},'MarkerEdgeColor',colour_line{ndataset})
            m{ndataset}.Color = colour_line{ndataset};
            %             hold on
            %             r(nshuffle) = rectangle('Position',pos)

            xlim([0 0.06])
            ylim([0 0.25])
            set(gca,'xscale','log')
            xticks(flip(p_val_threshold(low_threshold)))
            xticklabels({'0.001','0.002','0.005','0.01','0.02','0.05'})

            ylabel('Proportion of false events')
            hold on
            xlabel('Sequenceness p value')
            title(Behavioural_epoches{epoch})
            %             title('Multitrack event proportion vs p value')
        end
        ax = gca;
        ax.FontSize = 12;
        set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
        nfig = nfig +1;
    end
    legend([s(2),s(3),s(4),s(5)], {dataset_type{2},dataset_type{3},dataset_type{4},dataset_type{5}},'Location','northwest')
    sgtitle(shuffle_type{nmethod})

    cd ground_truth_original\Figure
    filename = sprintf('wcorr %s emprical false positive rates.pdf',shuffle_type{nmethod})
    saveas(gcf,filename)
    filename = sprintf('wcorr %s emprical false positive rates.fig',shuffle_type{nmethod})
    saveas(gcf,filename)
    cd ..
    cd ..
    close
end
%     legend([s(1),s(2),s(3),s(4)], {'spike train circular shift','place field circular shift','place bin circular shift','time bin permutation'},'Position',[0.35 0.2 0.05 0.05])
%     sgtitle('Multitrack event proportion vs Sequenceness p value')

end



