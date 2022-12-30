function [] = plot_ground_truth_multitracks_event_combined(folders,option,method)

total_number = 0;

for nfolder = 1:length(folders)

    cd(folders{nfolder})
    load decoded_replay_events
    total_number = total_number + length(decoded_replay_events(1).replay_events); % For cell id shuffles, get
    cd ..
end

index = [];
epoch_index = [];
for nmethod = 1:length(method)

    cd ground_truth_original
    if strcmp(method{nmethod},'wcorr')
        load  log_odd_wcorr_one_shuffle
        log_odd_compare{1} = log_odd
        load log_odd_wcorr_two_shuffles
        log_odd_compare{2} = log_odd
        load log_odd_wcorr_one_shuffle
        log_odd_compare{3} = log_odd
        load log_odd_wcorr
        log_odd_compare{4} = log_odd

    elseif strcmp(method{nmethod},'spearman')
        load log_odd_spearman_one_shuffle
        log_odd_compare{1} = log_odd
        load log_odd_spearman_two_shuffles
        log_odd_compare{2} = log_odd
        load log_odd_spearman_one_shuffle
        log_odd_compare{3} = log_odd
        load log_odd_spearman
        log_odd_compare{4} = log_odd

    elseif strcmp(method{nmethod},'linear')
        load log_odd_linear_one_shuffle
        log_odd_compare{1} = log_odd
        load log_odd_linear_two_shuffles
        log_odd_compare{2} = log_odd
        load log_odd_linear_one_shuffle
        log_odd_compare{3} = log_odd
        load log_odd_linear
        log_odd_compare{4} = log_odd

    elseif strcmp(method{nmethod},'path')
        load log_odd_path_one_shuffle
        log_odd_compare{1} = log_odd
        load log_odd_path_two_shuffles
        log_odd_compare{2} = log_odd
        load log_odd_path_one_shuffle
        log_odd_compare{3} = log_odd
        load log_odd_path
        log_odd_compare{4} = log_odd
    end
    cd ..

    % Get index for PRE,  RUN, POST
    % states = [-1 0 1 2 3 4 5]; % sleepPRE, awakePRE, RUN T1, RUN T2, sleepPOST and awakePOST and multi-track
    % states = [-1 0 1 2 3 4]; % sleepPRE, awakePRE, RUN T1, RUN T2, sleepPOST and awakePOST
    %     states = [-1 0 1 2];

    for nshuffle = 1:4
        if nshuffle == 3
            jump_index = intersect(find(log_odd_compare{nshuffle}.multi_track == 1),...;
                find(log_odd_compare{nshuffle}.max_jump <= 8));
            index{nmethod}{nshuffle}{1} = intersect(jump_index,find(log_odd_compare{nshuffle}.track==1));
            index{nmethod}{nshuffle}{2} = intersect(jump_index,find(log_odd_compare{nshuffle}.track==2));

        else
            index{nmethod}{nshuffle}{1} = intersect(find(log_odd_compare{nshuffle}.multi_track == 1),...
                find(log_odd_compare{nshuffle}.track==1));
            index{nmethod}{nshuffle}{2} = intersect(find(log_odd_compare{nshuffle}.multi_track == 1),...
                find(log_odd_compare{nshuffle}.track==2));
        end

    end




    % index = epoch_index;

    if strcmp(option,'common')
        for nshuffle = 1:4
            data{nmethod}{nshuffle} = log_odd_compare{nshuffle}.common_zscored.original(2,:);
            log_pval{nmethod}{nshuffle} = log_odd_compare{nshuffle}.pvalue;
        end
    elseif strcmp(option,'original')
        for nshuffle = 1:4
            data{nmethod}{nshuffle} = log_odd_compare{nshuffle}.normal_zscored.original(2,:);
            log_pval{nmethod}{nshuffle} = log_odd_compare{nshuffle}.pvalue;
        end
    elseif strcmp(option,'global remapped')
        data{nmethod} = log_odd_compare{nshuffle}.normal_zscored.global_remapped_original(2,:);
    end
end

log_odd_difference = [];
percent_sig_events = [];
for nmethod = 1:length(method)
    colour_line= {'b','r','g','k'};
    %     colour_line2 = {'k--','b--','r--','g--'};
    colour_symbol={'bo','ro','go','ko'};
    Behavioural_epoches = {'PRE','RUN','POST','Multi'};

    fig = figure(1)
    fig.Position = [834 116 850 700];

    %     text1 = sprintf('%s',method)
    filename = sprintf('%s multitracks shuffle comparisions',method{nmethod})
    sgtitle(filename);
    p_val_threshold = flip([-3:0.1:-1.3]);
    alpha_level = linspace(0.01,1,length(p_val_threshold));

    for nshuffle = 1:length(log_pval{nmethod})
        for threshold = 1:length(p_val_threshold)
            % find the event index with p value lower than the threshold
            %         current_index = intersect(find(log_pval(index{epoch})< p_val_threshold(threshold)),find(log_pval(index{epoch})> p_val_threshold(threshold+1)));
            %             current_index = find(log_pval{nmethod}{nshuffle}(index{nmethod}{nshuffle}{epoch})< p_val_threshold(threshold));
            track_1_index = find(log_pval{nmethod}{nshuffle}(index{nmethod}{nshuffle}{1})< p_val_threshold(threshold));
            track_2_index = find(log_pval{nmethod}{nshuffle}(index{nmethod}{nshuffle}{2})< p_val_threshold(threshold));

            log_odd_difference{nmethod}{nshuffle}(threshold) = mean(data{nmethod}{nshuffle}(index{nmethod}{nshuffle}{1}(track_1_index))) ...
                - mean(data{nmethod}{nshuffle}(index{nmethod}{nshuffle}{2}(track_2_index)));

            %             no_sig_events(threshold) = length(current_index)/length(data{nmethod});
            percent_sig_events{nmethod}{nshuffle}(threshold) = (length(track_1_index) + length(track_2_index)) /total_number;

            sc(nshuffle) = scatter(log_odd_difference{nmethod}{nshuffle}(threshold),percent_sig_events{nmethod}{nshuffle}(threshold),colour_line{nshuffle},'filled','MarkerFaceAlpha',alpha_level(threshold))
            hold on
            xlabel('Mean log odd difference')
            ylabel('Proportion of all replay events')
        end
        %         no_sig_events = [];
        %         mean_log_odd = [];
        %         SE_log_odd = [];
        ylim([0 0.4])

    end
    filename = sprintf('%s multitracks shuffle comparisions',method{nmethod})
    sgtitle(filename)
    legend([sc(1),sc(2),sc(3),sc(4)], {'1 shuffle','2 shuffles','1 shuffle + jump distance','3 shuffles'})

    cd ground_truth_original\Figure
    filename = sprintf('%s multitracks shuffle comparisions.pdf',method{nmethod})
    saveas(gcf,filename)
    clf
    cd ..
    cd ..
    % filename = sprintf('%s (%s) ground truth plot.pdf',method,option)
    % saveas(gcf,filename)
end


log_odd_difference = [];
log_odd_difference_CI = [];
percent_sig_events = [];


for nmethod = 1:length(method)
    filename = sprintf('%s multitracks shuffle comparisions',method{nmethod})
    sgtitle(filename);
    p_val_threshold = flip([-3:0.01:-1.3]); % fine resolution
    for nshuffle = 1:length(log_pval{nmethod})
        for threshold = 1:length(p_val_threshold)
            % find the event index with p value lower than the threshold
            %             current_index = find(log_pval{nmethod}(index{nmethod}{epoch})< p_val_threshold(threshold));
            difference = [];
            percent = [];

            parfor nboot = 1:1000 % Bootstrapping 1000 times
                track_1_index = datasample(find(log_pval{nmethod}{nshuffle}(index{nmethod}{nshuffle}{1})< p_val_threshold(threshold))...
                    ,length(find(log_pval{nmethod}{nshuffle}(index{nmethod}{nshuffle}{1})< p_val_threshold(threshold))));
                track_2_index = datasample(find(log_pval{nmethod}{nshuffle}(index{nmethod}{nshuffle}{2})< p_val_threshold(threshold))...
                    ,length(find(log_pval{nmethod}{nshuffle}(index{nmethod}{nshuffle}{2})< p_val_threshold(threshold))));

                % Original
                difference(nboot) = mean(data{nmethod}{nshuffle}(index{nmethod}{nshuffle}{1}(track_1_index))) ...
                    - mean(data{nmethod}{nshuffle}(index{nmethod}{nshuffle}{2}(track_2_index)));

                percent(nboot) = (length(track_1_index)+length(track_2_index))/total_number;
            end

            log_odd_difference{nmethod}{nshuffle}(threshold,:) = difference;
            log_odd_difference_CI{nmethod}{nshuffle}(threshold,:) = prctile(difference,[2.5 97.5]);
            percent_sig_events{nmethod}{nshuffle}(threshold) = percent(1);
        end

        x = mean(log_odd_difference{nmethod}{nshuffle},2)';
        y = percent_sig_events{nmethod}{nshuffle};
        y(isnan(x)) = [];
        x(isnan(x)) = [];
        UCI = log_odd_difference_CI{nmethod}{nshuffle}(:,2)';
        UCI(isnan(UCI)) = [];
        LCI = log_odd_difference_CI{nmethod}{nshuffle}(:,1)';
        LCI(isnan(LCI)) = [];

        if isempty(x) == 1
            continue
        else
            p(nshuffle) = patch([LCI fliplr(UCI)], [y fliplr(y)], colour_line{nshuffle},'FaceAlpha','0.25','LineStyle','none');
        end
        hold on
        plot(x,y,colour_line{nshuffle});

        xlabel('mean log odd difference')
        ylabel('Proportion of all replay events')
    end
    ylim([0 0.4])
    legend([p(1),p(2),p(3),p(4)], {'1 shuffle','2 shuffles','1 shuffle + jump distance','3 shuffles'})
    %             legend([sc(1),sc(2),sc(3),sc(4)], {'1 shuffle','2 shuffles','1 shuffle + jump distance','3 shuffles'})
    cd ground_truth_original\Figure
    filename = sprintf('%s (%s) shuffle multitracks comparisions CI.pdf',method{nmethod},option)
    saveas(gcf,filename)
    % filename = sprintf('%s (%s) ground truth plot.pdf',method,option)
    % saveas(gcf,filename)
    cd ..
    cd ..
    clf
    %
    %     fig = figure(2)
    %     fig.Position = [834 116 900 860];
    %     for epoch = 1:3
    %         subplot(2,2,epoch)
    %         %         hist_edges = ;
    %         % 3 shuffles at 0.05 and 2 shuffles
    %         [c index] = min(abs(percent_sig_events{nmethod}{2}{epoch}(:)-percent_sig_events{nmethod}{3}{epoch}(1))); % Get p value index at roughly equivalent propotion of events
    %         text1 = sprintf('3 shuffles p =< 0.05',10^p_val_threshold(1));
    %         text2 = sprintf('2 shuffles p =< %.3f',10^p_val_threshold(index));
    %
    %         h(1) = histogram(log_odd_difference{nmethod}{3}{epoch}(1,:),10,'FaceColor','green','FaceAlpha',0.3,'Normalization','probability');
    %         box off
    %         hold on
    %         h(2) = histogram(log_odd_difference{nmethod}{2}{epoch}(index,:),10,'FaceColor','red','FaceAlpha',0.3,'Normalization','probability');
    %         if log_odd_difference_CI{nmethod}{3}{epoch}(1,2) <= log_odd_difference_CI{nmethod}{2}{epoch}(index,1) || log_odd_difference_CI{nmethod}{3}{epoch}(1,1) >= log_odd_difference_CI{nmethod}{2}{epoch}(index,2)
    %             scatter(prctile(log_odd_difference{nmethod}{2}{epoch}(index,:),50),0.3,'r','*')
    %         end
    %         % 3 shuffles at 0.05 and 1 shuffle
    %         [c index] = min(abs(percent_sig_events{nmethod}{1}{epoch}(:)-percent_sig_events{nmethod}{3}{epoch}(1))); % Get p value index at roughly equivalent propotion of events
    %         text3 = sprintf('1 shuffles p =< %.3f',10^p_val_threshold(index));
    %         h(3) = histogram(log_odd_difference{nmethod}{1}{epoch}(index,:),10,'FaceColor','blue','FaceAlpha',0.3,'Normalization','probability');
    %         if log_odd_difference_CI{nmethod}{3}{epoch}(1,2) <= log_odd_difference_CI{nmethod}{1}{epoch}(index,1) || log_odd_difference_CI{nmethod}{3}{epoch}(1,1) >= log_odd_difference_CI{nmethod}{1}{epoch}(index,2)
    %             scatter(prctile(log_odd_difference{nmethod}{1}{epoch}(index,:),50),0.3,'b','*')
    %         end
    %
    %         xlabel('mean log odd difference')
    %         ylabel('Propotion')
    %         ylim([0 0.5])
    %
    %         titlename = sprintf('%s (propotion of events detected = %.3f)',Behavioural_epoches{epoch},percent_sig_events{nmethod}{3}{epoch}(1));
    %         title(titlename);
    %
    %         lgd = legend([h(1) h(2) h(3)],{text1,text2,text3},'Location','northwest');
    %     end
    %     titlename = sprintf('%s 3 shuffle vs other shuffles at equivalent proportion of replay events detected',method{nmethod});
    %     sgtitle(titlename)
    %     filename = sprintf('%s 3 shuffles multitracks comparisions histogram.pdf',method{nmethod});
    %     saveas(gcf,filename)
    %     clf
    %
    %     fig = figure(3)
    %     fig.Position = [834 116 900 860];
    %     for epoch = 1:3
    %         subplot(2,2,epoch)
    %         %         hist_edges = ;
    %         % 2 shuffles at 0.05 and 1 shuffle
    %         [c index] = min(abs(percent_sig_events{nmethod}{1}{epoch}(:)-percent_sig_events{nmethod}{2}{epoch}(1))); % Get p value index at roughly equivalent propotion of events
    %         text1 = sprintf('2 shuffles p =< 0.05',10^p_val_threshold(1));
    %         text2 = sprintf('1 shuffles p =< %.3f',10^p_val_threshold(index));
    %         h(1) = histogram(log_odd_difference{nmethod}{2}{epoch}(1,:),10,'FaceColor','red','FaceAlpha',0.3,'Normalization','probability');
    %         box off
    %         hold on
    %         h(2) = histogram(log_odd_difference{nmethod}{1}{epoch}(index,:),10,'FaceColor','blue','FaceAlpha',0.3,'Normalization','probability');
    %         if log_odd_difference_CI{nmethod}{2}{epoch}(1,2) <= log_odd_difference_CI{nmethod}{1}{epoch}(index,1) || log_odd_difference_CI{nmethod}{2}{epoch}(1,1) >= log_odd_difference_CI{nmethod}{1}{epoch}(index,2)
    %             scatter(prctile(log_odd_difference{nmethod}{1}{epoch}(index,:),50),0.3,'b','*')
    %         end
    %
    %         xlabel('mean log odd difference')
    %         ylabel('Propotion')
    %         ylim([0 0.5])
    %
    %         titlename = sprintf('%s (propotion of events detected = %.3f)',Behavioural_epoches{epoch},percent_sig_events{nmethod}{3}{epoch}(1));
    %         title(titlename);
    %
    %         lgd = legend([h(1) h(2)],{text1,text2},'Location','northwest');
    %     end
    %     titlename = sprintf('%s 2 shuffles vs 1 shuffle at equivalent proportion of replay events detected',method{nmethod});
    %     sgtitle(titlename)
    %     filename = sprintf('%s 2 shuffle multitracks comparisions histogram.pdf',method{nmethod});
    %     saveas(gcf,filename)
    %     clf

    cd ..
    cd ..
end

end

