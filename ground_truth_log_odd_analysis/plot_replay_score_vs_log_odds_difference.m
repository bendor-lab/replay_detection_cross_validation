function [] = plot_replay_score_vs_log_odds_difference(folders,option,method)

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
% load log_odd_difference_optimisation_original
% load log_odd_difference_optimisation
load log_odd_difference_optimisation_cross_experiment_shuffled
cd ..

nfig = 1;

%% Replay score vs mean log odds difference
p_val_threshold = round(10.^[-3:0.01:log10(0.2)],4);
% p_val_threshold = 0.001:0.0001:0.2;
p_val_threshold = flip(unique(p_val_threshold));
p_val_threshold(1) = 0.2;
low_threshold = find(ismembertol(p_val_threshold,[0.0501 0.02 0.01 0.005 0.002 0.001],eps) == 1);
method_type = {'wcorr 1 shuffle','wcorr 1 shuffle + jump distance','wcorr 2 shuffles','wcorr 3 shuffles'...
    'spearman median spike','spearman all spikes',...
    'linear 1 shuffle','linear 2 shuffles'};
% method_type = {'wcorr 2 shuffles',...
%     'spearman median spike','spearman all spikes',...
%     'linear 2 shuffles'};

% colour_line= {[34,94,168]/255,...
%     [253,141,60]/255,[254,217,11]/255,...
%     [227,26,28]/255};
colour_line= {[0,0,0]/255,[67,147,195]/255,[244,109,67]/255,[215,48,39]/255};
alpha_level = linspace(0.2,1,length(low_threshold));
size_level = linspace(40,80,length(low_threshold));

for nmethod = 1:length(method)
    fig = figure(nmethod)
    fig.Position = [834 116 850 700];
    for epoch = 1:3
        subplot(2,2,epoch)
        for nshuffle = 1:2

            y = mean(replay_score_compare{nmethod}{nshuffle}{epoch},2)';
            
            x = mean(log_odd_difference{nmethod}{nshuffle}{epoch},2)';
            UCI = log_odd_difference_CI{nmethod}{nshuffle}{epoch}(:,2)';
            LCI = log_odd_difference_CI{nmethod}{nshuffle}{epoch}(:,1)';
            
            % Because y does not always increase, need to re-sort the order
            [y id] = sort(y);
            x = x(id);
            UCI = UCI(id);
            LCI = LCI(id);

            if nshuffle == 1
                for threshold = 1:length(low_threshold)
                    scatter(x(find(id == low_threshold(threshold))),y(find(id == low_threshold(threshold))),size_level(threshold),'filled','MarkerFaceColor',colour_line{epoch},'MarkerFaceAlpha',alpha_level(threshold),'MarkerEdgeColor',colour_line{epoch})
                    hold on
                end
            else
                for threshold = 1:length(low_threshold)
                    scatter(x(find(id == low_threshold(threshold))),y(find(id == low_threshold(threshold))),size_level(threshold),'filled','MarkerFaceColor','k','MarkerFaceAlpha',alpha_level(threshold),'MarkerEdgeColor',colour_line{epoch})
                    hold on
                end
            end
            
            UCI(isnan(x)) = [];
            LCI(isnan(x)) = [];
            y(isnan(x)) = [];
            x(isnan(x)) = [];
            
            
            if nshuffle == 1
                p(nshuffle) = patch([LCI fliplr(UCI)], [y fliplr(y)], colour_line{epoch},'FaceAlpha','0.2','LineStyle','none');
                hold on
            else
                p(nshuffle) = patch([LCI fliplr(UCI)], [y fliplr(y)],'k','FaceAlpha','0.1','LineStyle','none');
                hold on
            end
            
%             if epoch == 2
                xlim([-0.5 4])
%             elseif epoch == 3
                xlim([-0.5 3])
%             end
            title(Behavioural_epoches{epoch});

            xlabel('mean log odds difference')
            ylabel('mean replay score')

        end
        %         legend([p(1),p(2),p(3)], {Behavioural_epoches{1},Behavioural_epoches{2},Behavioural_epoches{3}},'Position',[0.7 0.3 0.05 0.05])
        if nshuffle == 2
            legend([p(1),p(2)], {'original data','cell ID randomised'},'Location','southeast')
        else
            legend([p(1),p(3)], {'original data','cross experiment cell ID randomised'},'Location','southeast')
        end
        %         legend([p(1),p(3)], {'original data','shuffled data'},'Location','southeast')
        %         legend([p(1),p(2)], {'original data','shuffled data'},'Position',[0.7 0.3 0.05 0.05])
        set(gca,"TickDir","out",'box', 'off','Color','none')
        ax = gca;
        ax.FontSize = 14;
    end
    cd ground_truth_original\Figure
    if nshuffle == 2
        filename = sprintf('%s framework replay score vs log odds difference.pdf',method{nmethod})
        sgtitle(method{nmethod})
        saveas(gcf,filename)
    else
        filename = sprintf('%s framework replay score vs log odds difference cross experiment randomised.pdf',method{nmethod})
        sgtitle(method{nmethod})
        saveas(gcf,filename)
    end
    cd ..
    cd ..
    clf
end

% colour_line= {[34,94,168]/255,...
%     [253,141,60]/255,[254,217,11]/255,...
%     [227,26,28]/255};

alpha_level = linspace(0.2,1,length(low_threshold));
size_level = linspace(40,80,length(low_threshold));
colour_line= {[127,205,187]/255,[65,182,196]/255,[34,94,168]/255,[37,52,148]/255,...
    [253,141,60]/255,[254,217,11]/255,...
    [128,0,38]/255,[227,26,28]/255};

for nmethod = 1:length(method)
    fig = figure(nmethod)
    fig.Position = [834 116 850 700];
    for epoch = 1:3
        subplot(2,2,epoch)
        for nshuffle = [1 2]

            y = mean(replay_score_compare{nmethod}{nshuffle}{epoch},2)';
            x = mean(log_odd_difference{nmethod}{nshuffle}{epoch},2)';
            UCI = log_odd_difference_CI{nmethod}{nshuffle}{epoch}(:,2)';
            LCI = log_odd_difference_CI{nmethod}{nshuffle}{epoch}(:,1)';
            
            % Because y does not always increase, need to re-sort the order
            [y id] = sort(y);
            x = x(id);
            UCI = UCI(id);
            LCI = LCI(id);

            if nshuffle == 1
                for threshold = 1:length(low_threshold)
                    scatter(x(find(id == low_threshold(threshold))),y(find(id == low_threshold(threshold))),size_level(threshold),'filled','MarkerFaceColor',colour_line{nmethod},'MarkerFaceAlpha',alpha_level(threshold),'MarkerEdgeColor',colour_line{nmethod})
                    hold on
                end
            else
                for threshold = 1:length(low_threshold)
                    scatter(x(find(id == low_threshold(threshold))),y(find(id == low_threshold(threshold))),size_level(threshold),'filled','MarkerFaceColor','k','MarkerFaceAlpha',alpha_level(threshold),'MarkerEdgeColor',colour_line{nmethod})
                    hold on
                end
            end
            
            UCI(isnan(x)) = [];
            LCI(isnan(x)) = [];
            y(isnan(x)) = [];
            x(isnan(x)) = [];
            
            


            if nshuffle == 1
                p(nshuffle) = patch([LCI fliplr(UCI)], [y fliplr(y)], colour_line{nmethod},'FaceAlpha','0.2','LineStyle','none');
                hold on
            else
                p(nshuffle) = patch([LCI fliplr(UCI)], [y fliplr(y)],'k','FaceAlpha','0.1','LineStyle','none');
                hold on
            end

            title(Behavioural_epoches{epoch});

            xlabel('mean log odds difference')
            ylabel('mean replay score')

        end
        %         legend([p(1),p(2),p(3)], {Behavioural_epoches{1},Behavioural_epoches{2},Behavioural_epoches{3}},'Position',[0.7 0.3 0.05 0.05])
        legend([p(1),p(3)], {'original data','shuffled data'},'Location','southeast')
%         legend([p(1),p(2)], {'original data','shuffled data'},'Position',[0.7 0.3 0.05 0.05])
        set(gca,"TickDir","out",'box', 'off','Color','none')
        ax = gca;
        ax.FontSize = 14;
    end
    cd ground_truth_original\Figure
    if nshuffle == 2
        filename = sprintf('%s replay score vs log odds difference comparisions CI.pdf',method{nmethod})
    elseif nshuffle == 3
        filename = sprintf('%s replay score vs log odds difference comparisions CI cross experiment randmoised.pdf',method{nmethod})
    end
    sgtitle(method{nmethod})
    saveas(gcf,filename)
    cd ..
    cd ..
end

%% Demonstrate replay cross-validation framework

p_val_threshold = round(10.^[-3:0.01:log10(0.2)],4);
% p_val_threshold = 0.001:0.0001:0.2;
p_val_threshold = flip(unique(p_val_threshold));
p_val_threshold(1) = 0.2;
low_threshold = find(ismembertol(p_val_threshold,[0.0501 0.02 0.01 0.005 0.002 0.001],eps) == 1);
method_type = {'wcorr 1 shuffle','wcorr 1 shuffle + jump distance','wcorr 2 shuffles','wcorr 3 shuffles'...
    'spearman median spike','spearman all spikes',...
    'linear 1 shuffle','linear 2 shuffles'};
% method_type = {'wcorr 2 shuffles',...
%     'spearman median spike','spearman all spikes',...
%     'linear 2 shuffles'};

% colour_line= {[34,94,168]/255,...
%     [253,141,60]/255,[254,217,11]/255,...
%     [227,26,28]/255};
colour_line= {[0,0,0]/255,[67,147,195]/255,[244,109,67]/255,[215,48,39]/255};
alpha_level = linspace(0.2,1,length(low_threshold));
size_level = linspace(40,80,length(low_threshold))


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
    for nmethod = 3
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
for nmethod = 3
    for nshuffle = 2:3
        for epoch = 2:3
            y = mean(percent_shuffle_events{nmethod}{epoch}{nshuffle}(low_threshold,:),2);

     

            if nshuffle == 2
                s(epoch) = scatter(p_val_threshold(low_threshold),y,20,marker_shape{epoch},...
                    'MarkerFaceColor',colour_line{epoch},'MarkerEdgeColor','k','MarkerFaceAlpha','0.3')
                hold on
            elseif nshuffle == 3
                s(epoch) = scatter(p_val_threshold(low_threshold),y,20,marker_shape{epoch},...
                    'MarkerFaceColor',colour_line{epoch},'MarkerEdgeColor','r','MarkerFaceAlpha','0.3')
                hold on
            end
                            m{epoch} = errorbar(p_val_threshold(low_threshold),y,y-percent_shuffle_events_CI{nmethod}{epoch}{nshuffle}(low_threshold,1),...
                    percent_shuffle_events_CI{nmethod}{epoch}{nshuffle}(low_threshold,2)-y,'.','MarkerFaceColor',colour_line{epoch},'MarkerEdgeColor',colour_line{epoch})
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
end

% Sig proportion and log odd at p value 0.05
subplot(2,2,3)
for nmethod = 3
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
        text{epoch} = sprintf('%s p =< %.3f (proportion = %.3f)',Behavioural_epoches{epoch},p_val_threshold(low_threshold(1)),mean(percent_shuffle_events{nmethod}{epoch}{nshuffle}(low_threshold(1),:)));

        ylabel('Proportion of significant events')
        hold on
        xlabel('mean log odds difference')
        %             title(Behavioural_epoches{epoch})
        title('Original p value =< 0.05')
    end
end

% Sig proportion and log odd at mean false positive rate 0.05 （cell id shuffle）
%         tempt = find(mean(percent_multi_events{nmethod}{nshuffle}{epoch},2) <= 0.05);
for nmethod = 3
    for nshuffle = 2
        for epoch = 2:3
            [c index(epoch)] = min(abs(mean(percent_shuffle_events{nmethod}{epoch}{nshuffle},2) - 0.05));
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

            patch(error_x,error_y,colour_line{epoch},'EdgeColor','k','FaceAlpha','0.2','EdgeAlpha','0.7')
            textc{epoch} = sprintf('%s p =< %.3f (proportion = %.3f)',Behavioural_epoches{epoch},p_val_threshold(index(epoch)),mean(percent_shuffle_events{nmethod}{epoch}{nshuffle}(index(epoch),:)));

            ylabel('Proportion of significant events')
            hold on
            xlabel('mean log odds difference')
            %             title(Behavioural_epoches{epoch})
            title('Original vs cell id shuffle corrected p value')
        end
    end
end


% Sig proportion and log odd at mean false positive rate 0.05 （cross experiment shuffle）
for nmethod = 3
    for nshuffle = 3
        for epoch = 2:3
            [c index(epoch)] = min(abs(mean(percent_shuffle_events{nmethod}{epoch}{nshuffle},2) - 0.05));
            x = mean(mean(log_odd_difference{nmethod}{1}{epoch}(index(epoch),:),2),2);
            y = mean(percent_sig_events{nmethod}{1}{epoch}(index(epoch),:),2);
            sc3(epoch) = scatter(x,y,20,marker_shape{epoch},...
                'MarkerFaceColor',colour_line{epoch},'MarkerEdgeColor','r','MarkerFaceAlpha','0.5')
            %             hold on
            %             m{nshuffle} = errorbar(p_val_threshold(low_threshold),y,y-percent_multi_events_CI{nmethod}{nshuffle}{epoch}(low_threshold,1),...
            %                 percent_multi_events_CI{nmethod}{nshuffle}{epoch}(low_threshold,2)-y,'.','MarkerFaceColor',colour_line{nshuffle},'MarkerEdgeColor',colour_line{nshuffle})
            %             m{nshuffle}.Color = colour_line{nshuffle};
            hold on
            error_x = [log_odd_difference_CI{nmethod}{1}{epoch}(index(epoch),1) log_odd_difference_CI{nmethod}{1}{epoch}(index(epoch),2)...
                log_odd_difference_CI{nmethod}{1}{epoch}(index(epoch),2) log_odd_difference_CI{nmethod}{1}{epoch}(index(epoch),1)]
            error_y = [percent_sig_events_CI{nmethod}{1}{epoch}(index(epoch),1) percent_sig_events_CI{nmethod}{1}{epoch}(index(epoch),1)...
                percent_sig_events_CI{nmethod}{1}{epoch}(index(epoch),2) percent_sig_events_CI{nmethod}{1}{epoch}(index(epoch),2)]

            patch(error_x,error_y,colour_line{epoch},'EdgeColor','r','FaceAlpha','0.2','EdgeAlpha','0.7')
            textc3{epoch} = sprintf('%s p =< %.3f (proportion = %.3f)',Behavioural_epoches{epoch},p_val_threshold(index(epoch)),mean(percent_shuffle_events{nmethod}{epoch}{nshuffle}(index(epoch),:)));

            ylabel('Proportion of significant events')
            hold on
            xlabel('mean log odds difference')
            %             title(Behavioural_epoches{epoch})
            title('Original vs cell id vs cross experiment shuffles')
        end
    end
end
%         plot([0 0.1],[0 0.3],'r')
% legend([sc3(2),sc3(3)], {textc3{2},textc3{3}})

xlim([0.5 2.2])
ylim([0 0.45])
ax = gca;
ax.FontSize = 12;
set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
yline(0.1,'--','LineWidth',0.01);
% legend([s(2),s(3),sc(2),sc(3)], {text{2},text{3},textc{2},textc{3}})
legend([s(2),s(3),sc(2),sc(3),sc3(2),sc3(3)], {text{2},text{3},textc{2},textc{3},textc3{2},textc3{3}})
% legend([sc(2),sc(3),sc3(2),sc3(2)], {textc{2},textc{3},textc3{2},textc3{3}})

cd ground_truth_original\Figure
filename = 'log odds framework original vs cell id vs cross experiment randomised p value.pdf';
saveas(gcf,filename)
filename = 'log odds framework original vs cell id vs cross experiment randomised p value.fig';
saveas(gcf,filename)
cd ..
cd ..
clf
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

for nmethod = 1:2
    for epoch = 1:3
        fig = figure(nmethod+20)
        fig.Position = [834 116 850 700];
        subplot(2,2,epoch)

        for nshuffle = [1 3]
            x = mean(log_odd_difference{nmethod}{nshuffle}{epoch},2)';
            y = mean(percent_sig_events{nmethod}{nshuffle}{epoch},2)';

            if nshuffle == 1
                for threshold = 1:length(low_threshold)
                    scatter(x(low_threshold(threshold)),y(low_threshold(threshold)),'filled','MarkerFaceColor',colour_line{nmethod},'MarkerFaceAlpha',alpha_level(threshold),'MarkerEdgeColor',colour_line{nmethod})
                    hold on
                end
            elseif nshuffle >= 2
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
            elseif nshuffle >= 2
                p(nshuffle) = patch([LCI fliplr(UCI)], [y fliplr(y)],'k','FaceAlpha','0.2','LineStyle','none');
                hold on
            end
            xlabel('mean log odds difference')
            ylabel('Proportion of all replay events')
            title(Behavioural_epoches{epoch});
        end
        
    end

    ylim([0 1])

    if nshuffle == 2
        legend([p(1),p(2)], {'Original','Cell-id shuffle'},'Position',[0.7 0.3 0.05 0.05])
    elseif nshuffle == 3
        legend([p(1),p(3)], {'Original','cross experiment shuffle'},'Position',[0.7 0.3 0.05 0.05])
    end
    set(gca,"TickDir","out",'box', 'off','Color','none')
%     ax = gca;
%     ax.FontSize = 14;

%     cd ground_truth_original\Figure
%     if      nshuffle == 2
%         filename = sprintf('%s optimisation comparisions CI.pdf',method{nmethod});
%     elseif nshuffle == 3
%         filename = sprintf('%s optimisation comparisions CI cross experiment randomised.pdf',method{nmethod});
%     end
%     sgtitle(method{nmethod})
%     saveas(gcf,filename)
%     cd ..
%     cd ..
%     clf
end
%     legend([p(1) p(3) p(6)], {'Shuffle 3','Shuffle 2+3','Shuffle 1+2+3'})






end



