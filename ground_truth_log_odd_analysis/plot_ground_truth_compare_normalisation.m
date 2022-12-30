function plot_ground_truth_compare_normalisation(folders,option,method)
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

            elseif event_time >= time_range.track(1).behaviour(1) & event_time <= time_range.track(1).behaviour(2)
                total_number{nshuffle}(2) = total_number{nshuffle}(2) + 1*number_of_global_remapped_shuffles;
                session_total_number{nshuffle}(nfolder,2) = session_total_number{nshuffle}(nfolder,2) + 1*number_of_global_remapped_shuffles;

            elseif event_time >= time_range.track(2).behaviour(1) & event_time <= time_range.track(2).behaviour(2)
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
    if strcmp(method{nmethod},'wcorr normalised within')
        load log_odd_wcorr_PRE_place_POST_time_within
        log_odd_compare{nmethod}{1} = log_odd;
        load log_odd_wcorr_PRE_place_POST_time_within_global_remapped
        log_odd_compare{nmethod}{2} = log_odd;

    elseif strcmp(method{nmethod},'wcorr normalised across')
        load log_odd_wcorr_PRE_place_POST_time
        log_odd_compare{nmethod}{1} = log_odd;
        load log_odd_wcorr_PRE_place_POST_time_global_remapped
        log_odd_compare{nmethod}{2} = log_odd;

    elseif strcmp(method{nmethod},'linear normalised within')
        load log_odd_linear_PRE_place_POST_time_within
        log_odd_compare{nmethod}{1} = log_odd;
        load log_odd_linear_PRE_place_POST_time_within_global_remapped
        log_odd_compare{nmethod}{2} = log_odd;

    elseif strcmp(method{nmethod},'linear normalised across')
        load log_odd_linear_PRE_place_POST_time
        log_odd_compare{nmethod}{1} = log_odd;
        load log_odd_linear_PRE_place_POST_time_global_remapped
        log_odd_compare{nmethod}{2} = log_odd;

    end
    cd ..

    % Get index for PRE,  RUN, POST
    % states = [-1 0 1 2 3 4 5]; % sleepPRE, awakePRE, RUN T1, RUN T2, sleepPOST and awakePOST and multi-track
    % states = [-1 0 1 2 3 4]; % sleepPRE, awakePRE, RUN T1, RUN T2, sleepPOST and awakePOST
    states = [-1 0 1 2];
    
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
if exist('log_odd_difference_compare_normalisation.mat', 'file') ~= 2
    
    load log_odd_difference_optimisation
    % copy wcorr 2 shuffles
    log_odd_difference1{2} = log_odd_difference{3};
    log_odd_difference_CI1{2} = log_odd_difference_CI{3};
    percent_sig_events1{2} = percent_sig_events{3};
    percent_sig_events_CI1{2} = percent_sig_events_CI{3};
    percent_multi_events1{2} = percent_multi_events{3};
    percent_multi_events_CI1{2} = percent_multi_events_CI{3};
    percent_shuffle_events1{2} = percent_shuffle_events{3};
    percent_shuffle_events_CI1{2} = percent_shuffle_events_CI{3};

    % copy linear 2 shuffles
    log_odd_difference1{4} = log_odd_difference{8};
    log_odd_difference_CI1{4} = log_odd_difference_CI{8};
    percent_sig_events1{4} = percent_sig_events{8};
    percent_sig_events_CI1{4} = percent_sig_events_CI{8};
    percent_multi_events1{4} = percent_multi_events{8};
    percent_multi_events_CI1{4} = percent_multi_events_CI{8};
    percent_shuffle_events1{4} = percent_shuffle_events{8};
    percent_shuffle_events_CI1{4} = percent_shuffle_events_CI{8};
    
    clear log_odd_difference log_odd_difference_CI percent_sig_events percent_sig_events_CI ...
        percent_multi_events percent_multi_events_CI percent_shuffle_events percent_shuffle_events_CI

    for nmethod = [1 3]
        for epoch = 1:3
            for nshuffle = 1:length(log_pval{nmethod})
                tic
                
                parfor threshold = 1:length(p_val_threshold)
                    for nboot = 1:1000 % Bootstrapping 1000 times
                        s = RandStream('mrg32k3a','Seed',nboot);
                        
                        resampled_event = datasample(s,epoch_index{nmethod}{nshuffle}{epoch},length(epoch_index{nmethod}{nshuffle}{epoch}));
%                         resampled_event = datasample(epoch_index{nmethod}{nshuffle}{epoch},length(epoch_index{nmethod}{nshuffle}{epoch}));
                        
                        %     for track = 1:2
                        this_segment_event = resampled_event(find(segment_id{nmethod}{nshuffle}(1,resampled_event) == 1));
                        track_1_index = this_segment_event(find(log_pval{nmethod}{nshuffle}(1,this_segment_event) <= p_val_threshold(threshold)));
                        this_segment_event = resampled_event(find(segment_id{nmethod}{nshuffle}(1,resampled_event) > 1));
                        track_1_index = [track_1_index this_segment_event(find(log_pval{nmethod}{nshuffle}(1,this_segment_event) <= p_val_threshold(threshold)/2))];
                        
                        this_segment_event = resampled_event(find(segment_id{nmethod}{nshuffle}(2,resampled_event) == 1));
                        track_2_index = this_segment_event(find(log_pval{nmethod}{nshuffle}(2,this_segment_event) <= p_val_threshold(threshold)));
                        this_segment_event = resampled_event(find(segment_id{nmethod}{nshuffle}(2,resampled_event) > 1));
                        track_2_index = [track_2_index this_segment_event(find(log_pval{nmethod}{nshuffle}(2,this_segment_event) <= p_val_threshold(threshold)/2))];
                        
                        boot_replay_score(threshold,nboot) = mean([replay_score{nmethod}{nshuffle}(1,track_1_index) replay_score{nmethod}{nshuffle}(2,track_2_index)]);

                        multi_event_number = sum(ismember(track_1_index,track_2_index));
                        multi_event_percent = multi_event_number/(length(track_1_index)+length(track_2_index)-multi_event_number);
                        
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
    
    save log_odd_difference_compare_normalisation log_odd_difference log_odd_difference_CI percent_sig_events percent_sig_events_CI...
        percent_multi_events percent_multi_events_CI percent_shuffle_events percent_shuffle_events_CI replay_score_compare replay_score_compare_CI
    clear log_odd_difference1 log_odd_difference_CI1 percent_sig_events1 percent_sig_events_CI1 ...
        percent_multi_events1 percent_multi_events_CI1 percent_shuffle_events1 percent_shuffle_events_CI1
    clear boot_log_odd_difference boot_percent_sig_events boot_percent_multi_events boot_percent_shuffle_events
    cd ..
else
    load log_odd_difference_compare_normalisation
    cd ..
end
    
% 
% % colour_line = {[69,117,180]/255,[244,165,130]/255,[215,48,39]/255};
%  colour_line= {[8,81,156]/255,[67,147,195]/255,[244,109,67]/255,[215,48,39]/255};
%  p_val_threshold = round(10.^[-3:0.01:log10(0.2)],4);
% % p_val_threshold = 0.001:0.0001:0.2;
% p_val_threshold = flip(unique(p_val_threshold));
% p_val_threshold(1) = 0.2;
% p_value_to_plot = round(10.^[-3:0.1:log10(0.2)],4);
% p_value_to_plot = flip(unique(p_value_to_plot));
% p_value_to_plot(1) = 0.2;
% alpha_level = linspace(0.1,0.70,length(p_value_to_plot));
% low_threshold = find(ismembertol(p_val_threshold,p_value_to_plot,eps) == 1);
% 
% fig = figure(1)
% fig.Position = [834 116 850 700];
% for epoch = 1:3
%     subplot(2,2,epoch)
%     for nmethod = 1:4
%         x = log_odd_difference_original{nmethod}{1}{epoch}';
%         y = percent_sig_events_original{nmethod}{1}{epoch}';
% 
%         for threshold = 1:length(p_value_to_plot)
%             if p_value_to_plot(threshold) > 0.0501
%                 sc(nmethod) = scatter(x(low_threshold(threshold)),y(low_threshold(threshold)),'filled','MarkerFaceColor',colour_line{nmethod},'MarkerFaceAlpha',alpha_level(threshold))
%                 hold on
%             else
%                 sc(nmethod) = scatter(x(low_threshold(threshold)),y(low_threshold(threshold)),'filled','MarkerFaceColor',colour_line{nmethod},'MarkerFaceAlpha',alpha_level(threshold),'MarkerEdgeColor',colour_line{nmethod})
%                 hold on
%             end
%         end
%         %             plot(x,y,colour_line{nshuffle});
% 
% 
%         xlabel('mean log odd difference')
%         ylabel('Proportion of all replay events')
%     end
% 
%     ax = gca;
%     ax.FontSize = 12;
%     set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
%     title(Behavioural_epoches{epoch});
% 
%     ylim([0 1])
%     title(Behavioural_epoches{epoch});
% end
% legend([sc(1),sc(2),sc(3),sc(4)], {method{1},method{2},method{3},method{4}})
% 
% cd ground_truth_original\Figure
% filename = sprintf('normalisation comparision original.pdf')
% saveas(gcf,filename)
% clf
% cd ..
% cd ..


p_val_threshold = round(10.^[-3:0.01:log10(0.2)],4);
% p_val_threshold = 0.001:0.0001:0.2;
p_val_threshold = flip(unique(p_val_threshold));
p_val_threshold(1) = 0.2;
low_threshold = find(ismembertol(p_val_threshold,[0.0501 0.02 0.01 0.005 0.002 0.001],eps) == 1);
 colour_line= {[8,81,156]/255,[67,147,195]/255,[244,109,67]/255,[215,48,39]/255};

alpha_level = linspace(0.2,0.7,length(low_threshold));
fig = figure(2)
fig.Position = [834 116 850 700];
for epoch = 1:3
    subplot(2,2,epoch)
    for nmethod = 1:4
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
    ax = gca;
    ax.FontSize = 12;
    set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
    ylim([0 0.62])
    title(Behavioural_epoches{epoch});
end
%     legend([p(1) p(3) p(6)], {'Shuffle 3','Shuffle 2+3','Shuffle 1+2+3'})
legend([p(1),p(2),p(3),p(4)], {method{1},method{2},method{3},method{4}},'Position',[0.7 0.3 0.05 0.05])

cd ground_truth_original\Figure
filename = sprintf('normalisation comparisions CI.pdf')
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

fig = figure(2)
fig.Position = [834 116 850 885];
for nmethod = 1:4
    for epoch = 1:3
        subplot(2,2,epoch)
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
        ylim([0 0.1])
        ylabel('Proportion of cell-id shuffled events')
        hold on
        xlabel('Sequenceness p value')
        title(Behavioural_epoches{epoch})
        %             title('Multitrack event proportion vs p value')
        ax = gca;
        ax.FontSize = 12;
        set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
    end

end
legend([s(1),s(2),s(3),s(4)], {method{1},method{2},method{3},method{4}},'Position',[0.7 0.3 0.05 0.05])

%     legend([s(1),s(2),s(3),s(4)], {'spike train circular shift','place field circular shift','place bin circular shift','time bin permutation'},'Position',[0.35 0.2 0.05 0.05])
%     sgtitle('Multitrack event proportion vs Sequenceness p value')
cd ground_truth_original\Figure
filename = sprintf('normalisation comparision cell-id shuffled sig proportion.pdf')
saveas(gcf,filename)
cd ..
cd ..
clf



% Sig proportion and log odd at p value 0.05
fig = figure(3)
fig.Position = [834 116 850 885];
subplot(2,2,1)
for nmethod = 1:4
    for epoch = 1:3
        x = mean(mean(log_odd_difference{nmethod}{1}{epoch}(low_threshold(1),:),2),2);
        y = mean(percent_sig_events{nmethod}{1}{epoch}(low_threshold(1),:),2);
        s(nmethod) = scatter(x,y,20,marker_shape{epoch},...
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
        text{nmethod} = sprintf('%s p =< %.3f (proportion = %.3f)',method{nmethod},p_val_threshold(low_threshold(1)),mean(percent_shuffle_events{nmethod}{epoch}(low_threshold(1),:)));

        ylabel('Proportion of significant events')
        hold on
        xlabel('mean log odd difference')
        %             title(Behavioural_epoches{epoch})
        title('Original p value =< 0.05')
    end
end
legend([s(1),s(2),s(3),s(4)], {text{1},text{2},text{3},text{4}})

xlim([-1 3])
ylim([0 0.7])

% Sig proportion and log odd at shuffle-corrected p value 0.05
subplot(2,2,2)
%         tempt = find(mean(percent_multi_events{nmethod}{nshuffle}{epoch},2) <= 0.05);
for nmethod = 1:4
    for epoch = 1:3
        [c index(nmethod)] = min(abs(mean(percent_shuffle_events{nmethod}{epoch},2) - 0.05));
%          [c index(nmethod)] = min(abs(mean(percent_multi_events{nmethod}{1}{epoch},2) - 0.05));
        x = mean(mean(log_odd_difference{nmethod}{1}{epoch}(index(nmethod),:),2),2);
        y = mean(percent_sig_events{nmethod}{1}{epoch}(index(nmethod),:),2);
        s(nmethod) = scatter(x,y,20,marker_shape{epoch},...
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
        text{nmethod} = sprintf('%s p =< %.3f (proportion = %.3f)',method{nmethod},p_val_threshold(index(nmethod)),mean(percent_shuffle_events{nmethod}{epoch}(index(nmethod),:)));
%         text{nmethod} = sprintf('%s p =< %.3f (proportion = %.3f)',shuffle_type{nmethod},p_val_threshold(index(nmethod)),mean(percent_multi_events{nmethod}{1}{epoch}(index(nmethod),:)));

        ylabel('Proportion of significant events')
        hold on
        xlabel('mean log odd difference')
        %             title(Behavioural_epoches{epoch})
        title('Equivalent p value when cell-id shuffled sig proportion = 0.05')
    end
end
%         plot([0 0.1],[0 0.3],'r')
legend([s(1),s(2),s(3),s(4)], {text{1},text{2},text{3},text{4}})
xlim([-1 3])
ylim([0 0.7])

cd ground_truth_original\Figure
filename = sprintf('normalisation original vs shuffle-corrected p value.pdf')
saveas(gcf,filename)
cd ..
cd ..
clf



p_val_threshold = round(10.^[-3:0.01:log10(0.2)],4);
% p_val_threshold = 0.001:0.0001:0.2;
p_val_threshold = flip(unique(p_val_threshold));
p_val_threshold(1) = 0.2;
low_threshold = find(ismembertol(p_val_threshold,[0.0501 0.02 0.01 0.005 0.002 0.001],eps) == 1);

alpha_level = linspace(0.2,0.7,length(low_threshold));

for nmethod = [1 3]
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
            xlabel('mean log odd difference')
            ylabel('Proportion of all replay events')
            title(Behavioural_epoches{epoch});
        end
    end

    ylim([0 1])
    
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




p_val_threshold = round(10.^[-3:0.01:log10(0.2)],4);
p_val_threshold = flip(unique(p_val_threshold));
p_val_threshold(1) = 0.2;
low_threshold = find(ismembertol(p_val_threshold,[0.0501 0.02 0.01 0.005 0.002 0.001],eps) == 1);
Behavioural_epoches = {'PRE','RUN','POST'};
marker_shape = {'o','^','s'};

fig = figure(1)
fig.Position = [834 116 850 885];
subplot(2,2,1)
for nmethod = 1:length(method)
    for epoch = 1:3
        % shuffle-subtracted mean log odd difference
        shuffle_subtracted_log_odd_difference_CI{nmethod}{epoch} = prctile(log_odd_difference{nmethod}{1}{epoch} - log_odd_difference{nmethod}{2}{epoch},[2.5 97.5],2);
        shuffle_subtracted_log_odd_difference{nmethod}{epoch} = mean(log_odd_difference{nmethod}{1}{epoch} - log_odd_difference{nmethod}{2}{epoch},2);

        x_CI = shuffle_subtracted_log_odd_difference_CI{nmethod}{epoch}(low_threshold(1),:);
        x = shuffle_subtracted_log_odd_difference{nmethod}{epoch}(low_threshold(1),:);

        y = mean(percent_sig_events{nmethod}{1}{epoch}(low_threshold(1),:),2);
        s(nmethod) = scatter(x,y,20,marker_shape{epoch},...
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
        text{nmethod} = sprintf('%s p =< %.3f (proportion = %.3f)',method{nmethod},p_val_threshold(low_threshold(1)),mean(percent_shuffle_events{nmethod}{epoch}(low_threshold(1),:)));

        ylabel('Proportion of significant events')
        hold on
        xlabel('shuffle-subtracted mean log odd difference')
        %             title(Behavioural_epoches{epoch})
        title('Original p value =< 0.05')
    end
end
legend([s(1),s(2),s(3),s(4)],...
    {text{1},text{2},text{3},text{4}})

xlim([-1 3])
ylim([0 0.7])

% Sig proportion and log odd at shuffle-corrected p value 0.05
subplot(2,2,2)
%         tempt = find(mean(percent_multi_events{nmethod}{nshuffle}{epoch},2) <= 0.05);
for nmethod = 1:length(method)
    for epoch = 1:3
        [c index(nmethod)] = min(abs(mean(percent_shuffle_events{nmethod}{epoch},2) - 0.05));
%          [c index(nmethod)] = min(abs(mean(percent_multi_events{nmethod}{1}{epoch},2) - 0.05));
        x_CI = shuffle_subtracted_log_odd_difference_CI{nmethod}{epoch}(index(nmethod),:);
        x = shuffle_subtracted_log_odd_difference{nmethod}{epoch}(index(nmethod),:); 

        y = mean(percent_sig_events{nmethod}{1}{epoch}(index(nmethod),:),2);
        s(nmethod) = scatter(x,y,20,marker_shape{epoch},...
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
        text{nmethod} = sprintf('%s p =< %.3f (proportion = %.3f)',method{nmethod},p_val_threshold(index(nmethod)),mean(percent_shuffle_events{nmethod}{epoch}(index(nmethod),:)));
%         text{nmethod} = sprintf('%s p =< %.3f (proportion = %.3f)',shuffle_type{nmethod},p_val_threshold(index(nmethod)),mean(percent_multi_events{nmethod}{1}{epoch}(index(nmethod),:)));

        ylabel('Proportion of significant events')
        hold on
        xlabel('shuffle-subtracted mean log odd difference')
        %             title(Behavioural_epoches{epoch})
        title('Equivalent p value when cell-id shuffled sig proportion = 0.05')
    end
end
%         plot([0 0.1],[0 0.3],'r')
legend([s(1),s(2),s(3),s(4)],...
    {text{1},text{2},text{3},text{4}})

xlim([-1 3])
ylim([0 0.75])

cd ground_truth_original\Figure
filename = sprintf('normalisation comparision (shuffle-corrected) original vs shuffle-corrected p value.pdf')
saveas(gcf,filename)
filename = sprintf('normalisation comparision (shuffle-corrected) original vs shuffle-corrected p value.fig')
saveas(gcf,filename)
cd ..
cd ..
clf


end





