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
% load log_odd_difference_optimisation_original
load log_odd_difference_optimisation
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

for nmethod = 3
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
            elseif nshuffle == 2
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
            elseif nshuffle == 2
                p(nshuffle) = patch([LCI fliplr(UCI)], [y fliplr(y)],'k','FaceAlpha','0.1','LineStyle','none');
                hold on
            end

            title(Behavioural_epoches{epoch});

            xlabel('mean log odds difference')
            ylabel('mean replay score')

        end
        %         legend([p(1),p(2),p(3)], {Behavioural_epoches{1},Behavioural_epoches{2},Behavioural_epoches{3}},'Position',[0.7 0.3 0.05 0.05])
        legend([p(1),p(2)], {'original data','shuffled data'},'Location','southeast')
%         legend([p(1),p(2)], {'original data','shuffled data'},'Position',[0.7 0.3 0.05 0.05])

    end
    cd ground_truth_original\Figure
    filename = sprintf('%s framework replay score vs log odds difference.pdf',method{nmethod})
    sgtitle(method{nmethod})
    saveas(gcf,filename)
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
                    scatter(x(find(id == low_threshold(threshold))),y(find(id == low_threshold(threshold))),size_level(threshold),'filled','MarkerFaceColor',colour_line{nmethod},'MarkerFaceAlpha',alpha_level(threshold),'MarkerEdgeColor',colour_line{nmethod})
                    hold on
                end
            elseif nshuffle == 2
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
            elseif nshuffle == 2
                p(nshuffle) = patch([LCI fliplr(UCI)], [y fliplr(y)],'k','FaceAlpha','0.1','LineStyle','none');
                hold on
            end

            title(Behavioural_epoches{epoch});

            xlabel('mean log odds difference')
            ylabel('mean replay score')

        end
        %         legend([p(1),p(2),p(3)], {Behavioural_epoches{1},Behavioural_epoches{2},Behavioural_epoches{3}},'Position',[0.7 0.3 0.05 0.05])
        legend([p(1),p(2)], {'original data','shuffled data'},'Location','southeast')
%         legend([p(1),p(2)], {'original data','shuffled data'},'Position',[0.7 0.3 0.05 0.05])

    end
    cd ground_truth_original\Figure
    filename = sprintf('%s replay score vs log odds difference comparisions CI.pdf',method{nmethod})
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






end



