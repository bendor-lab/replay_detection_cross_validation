function [] = plot_ground_truth_compare_global_remapped_combined(folders,method,option)

total_number{1}(1:3) = 0;
total_number{2}(1:3) = 0;

for nfolder = 1:length(folders) % Get total number of replay events
    cd(folders{nfolder})

    load decoded_replay_events

    [~,time_range]=sort_replay_events([],'wcorr');

    for event = 1:length(decoded_replay_events(1).replay_events)

        event_time = decoded_replay_events(1).replay_events(event).timebins_edges(1);

        if event_time <= time_range.pre(2) %If PRE

            total_number{1}(1) = total_number{1}(1) + 1;

        elseif event_time >= time_range.post(1) %If POST

            total_number{1}(3) = total_number{1}(3) + 1;

        elseif event_time >= time_range.track(1).behaviour(1) & event_time <= time_range.track(1).behaviour(2)
            total_number{1}(2) = total_number{1}(2) + 1;
        elseif event_time >= time_range.track(2).behaviour(1) & event_time <= time_range.track(2).behaviour(2)
            total_number{1}(2) = total_number{1}(2) + 1;


        end
    end
    
    total_number{2}(1:3) = total_number{2}(1:3) + length(decoded_replay_events(1).replay_events); % For cell id shuffles, get 
    cd ..
end


index = [];
states = [-1 0 1 2];
epoch_index = [];

for nmethod = 1:length(method)
    log_odd_compare = [];
    cd ground_truth_original
    if strcmp(method{nmethod},'wcorr')

        load log_odd_wcorr
        log_odd_compare{1} = log_odd

        if exist('log_odd_wcorr_global_remapped.mat')~=2
            cd ..
            [log_odd] = extract_ground_truth_info(folders,'global remapped',method{nmethod},3,3)
            cd ground_truth_original
            save log_odd_wcorr_global_remapped log_odd

        else
            load log_odd_wcorr_global_remapped
        end

    elseif strcmp(method{nmethod},'spearman')
        load log_odd_spearman
        log_odd_compare{1} = log_odd

        if exist('log_odd_spearman_global_remapped.mat')~=2
            cd ..
            [log_odd] = extract_ground_truth_info(folders,'global remapped',method{nmethod},3,3)
            cd ground_truth_original
            save log_odd_spearman_global_remapped log_odd

        else
            load log_odd_spearman_global_remapped
        end

    elseif strcmp(method{nmethod},'linear')
        load log_odd_linear
        log_odd_compare{1} = log_odd

        if exist('log_odd_linear_global_remapped.mat')~=2
            [log_odd] = extract_ground_truth_info(folders,'global remapped',method{nmethod},3,3)
            save log_odd_linear_global_remapped log_odd
        else
            load log_odd_linear_global_remapped

        end

    elseif strcmp(method{nmethod},'path')
        load log_odd_path
        log_odd_compare{1} = log_odd

        if exist('log_odd_path_global_remapped.mat')~=2
            [log_odd] = extract_ground_truth_info(folders,'global remapped',method{nmethod},3,3)
            save log_odd_path_global_remapped log_odd
        else
            load log_odd_path_global_remapped
        end

    end
    log_odd_compare{2} = log_odd

    cd ..

    % Get index for PRE,  RUN, POST
    % states = [-1 0 1 2 3 4 5]; % sleepPRE, awakePRE, RUN T1, RUN T2, sleepPOST and awakePOST and multi-track
    % states = [-1 0 1 2 3 4]; % sleepPRE, awakePRE, RUN T1, RUN T2, sleepPOST and awakePOST


    for nshuffle = 1:2
        if nshuffle == 1 % Original
            for k=1:length(states) % sort track id and behavioral states (pooled from 5 sessions)
                state_index = find(log_odd_compare{nshuffle}.behavioural_state==states(k));

                if k == 1 % PRE
                    epoch_index{nmethod}{nshuffle}{1}{1} = intersect(state_index,find(log_odd_compare{nshuffle}.track==1));
                    epoch_index{nmethod}{nshuffle}{1}{2} = intersect(state_index,find(log_odd_compare{nshuffle}.track==2));
                elseif k == 2 % POST
                    epoch_index{nmethod}{nshuffle}{3}{1} = intersect(state_index,find(log_odd_compare{nshuffle}.track==1));
                    epoch_index{nmethod}{nshuffle}{3}{2} = intersect(state_index,find(log_odd_compare{nshuffle}.track==2));
                elseif k == 3 % RUN Track 1
                    epoch_index{nmethod}{nshuffle}{2}{1} = intersect(state_index,find(log_odd_compare{nshuffle}.track==1)); % Track 1 replay on Track 1 considered RUN Replay
                elseif k == 4 % RUN Track 2
                    epoch_index{nmethod}{nshuffle}{2}{2} = intersect(state_index,find(log_odd_compare{nshuffle}.track==2)); % Track 2 replay on Track 2 considered RUN Replay
                end
            end

            for epoch = 1:3
                %         tempt1 = datasample(epoch_index{epoch}{1},1000);
                %         tempt2 = datasample(epoch_index{epoch}{2},1000);

                index{nmethod}{nshuffle}{epoch} = [epoch_index{nmethod}{nshuffle}{epoch}{1} epoch_index{nmethod}{nshuffle}{epoch}{2}]
                %         index{nmethod}{epoch} = [ tempt1 tempt2 ]
            end
        elseif nshuffle == 2 % If cell id shuffled, behavioural epoches not important

            epoch_index{nmethod}{nshuffle}{1}{1} = find(log_odd_compare{nshuffle}.track==1);
            epoch_index{nmethod}{nshuffle}{1}{2} = find(log_odd_compare{nshuffle}.track==2);

            epoch_index{nmethod}{nshuffle}{2}{1} = find(log_odd_compare{nshuffle}.track==1);
            epoch_index{nmethod}{nshuffle}{2}{2} = find(log_odd_compare{nshuffle}.track==2);

            epoch_index{nmethod}{nshuffle}{3}{1} = find(log_odd_compare{nshuffle}.track==1);
            epoch_index{nmethod}{nshuffle}{3}{2} = find(log_odd_compare{nshuffle}.track==2);

            for epoch = 1:3

                index{nmethod}{nshuffle}{epoch} = [epoch_index{nmethod}{nshuffle}{epoch}{1} epoch_index{nmethod}{nshuffle}{epoch}{2}]
                %         index{nmethod}{epoch} = [ tempt1 tempt2 ]
            end

        end
    end
    % index = epoch_index;


    if strcmp(option,'original')
        data{nmethod}{1} = log_odd_compare{1}.normal_zscored.original(2,:); % Where (2,:) is T1/T2 log odd and (1,:) is T/F log odd
        data{nmethod}{2} = log_odd_compare{2}.normal_zscored.global_remapped_original(2,:);
    elseif strcmp(option,'common')
        data{nmethod}{1} = log_odd_compare{1}.common_zscored.original(2,:);
        data{nmethod}{2} = log_odd_compare{2}.common_zscored.global_remapped_original(2,:);
    end

    for nshuffle = 1:2
        log_pval{nmethod}{nshuffle} = log_odd_compare{nshuffle}.pvalue;
    end


    colour_line= {'b','r','g','k'};
    %     colour_line2 = {'k--','b--','r--','g--'};
    colour_symbol={'bo','ro','go','ko'};
    Behavioural_epoches = {'PRE','RUN','POST','Multi'};


    fig = figure(1)
    fig.Position = [834 116 850 700];

    %     text1 = sprintf('%s',method)
    filename = sprintf('%s (%s) Orignal vs shuffle.pdf',method{nmethod},option)
    sgtitle(filename);
    p_val_threshold = flip([-3:0.1:-1.3]);
    alpha_level = linspace(0.01,1,length(p_val_threshold));

    for epoch = 1:3
        subplot(2,2,epoch)

        for nshuffle = 1:2
  
            for threshold = 1:length(p_val_threshold)
                % find the event index with p value lower than the threshold
%                 current_index = find(log_pval{nmethod}{nshuffle}(index{nmethod}{nshuffle}{epoch})< p_val_threshold(threshold));
                track_1_index = find(log_pval{nmethod}{nshuffle}(epoch_index{nmethod}{nshuffle}{epoch}{1})< p_val_threshold(threshold));
                track_2_index = find(log_pval{nmethod}{nshuffle}(epoch_index{nmethod}{nshuffle}{epoch}{2})< p_val_threshold(threshold));

                % Original
                log_odd_difference{nmethod}{epoch}{nshuffle}(threshold) = mean(data{nmethod}{nshuffle}(epoch_index{nmethod}{nshuffle}{epoch}{1}(track_1_index))) ...
                    - mean(data{nmethod}{nshuffle}(epoch_index{nmethod}{nshuffle}{epoch}{2}(track_2_index)));

                percent_sig_events{nmethod}{epoch}{nshuffle}(threshold) = (length(track_1_index) + length(track_2_index))/total_number{nshuffle}(epoch);
                %             sc(nshuffle) = scatter(mean_log_odd{nmethod}{nshuffle}(threshold),percent_sig_events{nmethod}{nshuffle}(threshold),colour_line{nshuffle},'filled','MarkerFaceAlpha',alpha_level(threshold))
                %         plot(mean_log_odd,no_sig_events,colour_line{epoch})
                % Shuffle


                %             hold on
            end


            %     lgd.FontSize = 12;


            %     legend([sc(1),sc(2),sc(3),sc(4)], {'wcorr','spearman','linear','path'})
            %         title(Behavioural_epoches{epoch});

        end

    end


    for epoch = 1:3
        subplot(2,2,epoch)
        for threshold = 1:length(p_val_threshold)
            for n = 1:2
                %                 if no_of_events{nshuffle}{epoch}(threshold) >= 10
                if n == 1
                    sc(n) = scatter(log_odd_difference{nmethod}{epoch}{1}(threshold),percent_sig_events{nmethod}{epoch}{1}(threshold),colour_line{n},'filled','MarkerFaceAlpha',alpha_level(threshold))

                elseif n == 2
                    sc(n) = scatter(log_odd_difference{nmethod}{epoch}{2}(threshold),percent_sig_events{nmethod}{epoch}{2}(threshold),colour_line{n},'filled','MarkerFaceAlpha',alpha_level(threshold))
                end

                %                 end

                hold on
                xlabel('mean log odd difference')
                ylabel('Proportion of all replay events')
            end
        end
        ylim([0 0.4])
        title(Behavioural_epoches{epoch});
        if epoch ==1
            legend([sc(1),sc(2)], {'Original','Cell ID shuffle'})
        end

    end

    cd ground_truth_original\Figure

    filename = sprintf('%s (%s) Orignal vs shuffle combined.pdf',method{nmethod},option)
    saveas(gcf,filename)
    % filename = sprintf('%s (%s) ground truth plot.pdf',method,option)
    % saveas(gcf,filename)
    cd ..
    cd ..
    clf(fig)
end

p_val_threshold = flip([-3:0.01:-1.3]); % fine resolution
log_odd_difference = [];
log_odd_difference_CI = [];
percent_sig_events = [];

for nmethod = 1:length(method)
    for epoch = 1:3
        subplot(2,2,epoch)
        for nshuffle = 1:length(log_pval{nmethod})
            for threshold = 1:length(p_val_threshold)
                % find the event index with p value lower than the threshold
                %             current_index = find(log_pval{nmethod}(index{nmethod}{epoch})< p_val_threshold(threshold));
                difference = [];
                percent = [];

                parfor nboot = 1:10000 % Bootstrapping 1000 times
                    track_1_index = datasample(find(log_pval{nmethod}{nshuffle}(epoch_index{nmethod}{nshuffle}{epoch}{1})< p_val_threshold(threshold))...
                        ,length(find(log_pval{nmethod}{nshuffle}(epoch_index{nmethod}{nshuffle}{epoch}{1})< p_val_threshold(threshold))));
                    track_2_index = datasample(find(log_pval{nmethod}{nshuffle}(epoch_index{nmethod}{nshuffle}{epoch}{2})< p_val_threshold(threshold))...
                        ,length(find(log_pval{nmethod}{nshuffle}(epoch_index{nmethod}{nshuffle}{epoch}{2})< p_val_threshold(threshold))));

                    % Original
                    difference(nboot) = mean(data{nmethod}{nshuffle}(epoch_index{nmethod}{nshuffle}{epoch}{1}(track_1_index))) ...
                        - mean(data{nmethod}{nshuffle}(epoch_index{nmethod}{nshuffle}{epoch}{2}(track_2_index)));

                    percent(nboot) = (length(track_1_index)+length(track_2_index))/total_number{nshuffle}(epoch);
                end

                log_odd_difference{nmethod}{nshuffle}{epoch}(threshold,:) = difference;
                log_odd_difference_CI{nmethod}{nshuffle}{epoch}(threshold,:) = prctile(difference,[2.5 97.5]);
                percent_sig_events{nmethod}{nshuffle}{epoch}(threshold) = percent(1);
            end

            x = mean(log_odd_difference{nmethod}{nshuffle}{epoch},2)';
            y = percent_sig_events{nmethod}{nshuffle}{epoch};
            y(isnan(x)) = [];
            x(isnan(x)) = [];
            UCI = log_odd_difference_CI{nmethod}{nshuffle}{epoch}(:,2)';
            UCI(isnan(UCI)) = [];
            LCI = log_odd_difference_CI{nmethod}{nshuffle}{epoch}(:,1)';
            LCI(isnan(LCI)) = [];
            p(nshuffle) = patch([LCI fliplr(UCI)], [y fliplr(y)], colour_line{nshuffle},'FaceAlpha','0.25','LineStyle','none')
            hold on
            plot(x,y,colour_line{nshuffle});

            xlabel('mean log odd difference')
            ylabel('Proportion of all replay events')
        end

        ylim([0 0.4])
        title(Behavioural_epoches{epoch});
        if epoch ==1
            legend([p(1),p(2)], {'Original','Cell ID shuffle'})
        end
    end
    %
    %     for epoch = 1:3
    %         subplot(2,2,epoch)
    %         for threshold = 1:length(p_val_threshold)
    %             for n = 1:2
    %                 %                 if no_of_events{nshuffle}{epoch}(threshold) >= 10
    %                 if n == 1
    %                     sc(n) = scatter(log_odd_difference{nmethod}{epoch}{1}(threshold),percent_sig_events{nmethod}{epoch}{1}(threshold),colour_line{n},'filled','MarkerFaceAlpha',alpha_level(threshold))
    %
    %                 elseif n == 2
    %                     sc(n) = scatter(log_odd_difference{nmethod}{epoch}{2}(threshold),percent_sig_events{nmethod}{epoch}{2}(threshold),colour_line{n},'filled','MarkerFaceAlpha',alpha_level(threshold))
    %                 end
    %
    %                 %                 end
    %
    %                 hold on
    %                 xlabel('mean log odd difference')
    %                 ylabel('Proportion of all replay events')
    %             end
    %         end
    %         ylim([0 0.4])
    %         title(Behavioural_epoches{epoch});
    %         if epoch ==1
    %             legend([sc(1),sc(2)], {'Original','Cell ID shuffle'})
    %         end
    %
    %     end
    %
    cd ground_truth_original\Figure

    filename = sprintf('%s (%s) Orignal vs shuffle combined CI.pdf',method{nmethod},option)
    saveas(gcf,filename)
    % filename = sprintf('%s (%s) ground truth plot.pdf',method,option)
    % saveas(gcf,filename)
    cd ..
    cd ..
    clf(fig)

    %     fig = figure(2)
    %     fig.Position = [834 116 900 860];
    %     for epoch = 1:3
    %         subplot(2,2,epoch)
    %         %         hist_edges = ;
    %
    %         % Spearman at 0.05 vs Wcorr
    %         [c index] = min(abs(percent_sig_events{1}{epoch}(:)-percent_sig_events{2}{epoch}(1))); % Get p value index at roughly equivalent propotion of events
    %
    %         h(1) = histogram(log_odd_difference{2}{epoch}(1,:),10,'FaceColor','red','FaceAlpha',0.3,'Normalization','probability');
    %         hold on
    %         h(2) = histogram(log_odd_difference{1}{epoch}(index,:),10,'FaceColor','blue','FaceAlpha',0.3,'Normalization','probability');
    %         xlabel('mean log odd difference')
    %         ylabel('Propotion')
    %         ylim([0 0.4])
    %         if log_odd_difference_CI{2}{epoch}(1,2) <= log_odd_difference_CI{1}{epoch}(index,1) || log_odd_difference_CI{2}{epoch}(1,1) >= log_odd_difference_CI{1}{epoch}(index,2)
    %             scatter(prctile(log_odd_difference{1}{epoch}(index,:),50),0.3,'b','*')
    %         end
    %
    %         titlename = sprintf('%s (propotion of events detected = %.3f)',Behavioural_epoches{epoch},percent_sig_events{2}{epoch}(1));
    %         title(titlename);
    %
    %         text1 = 'Spearman p =< 0.05';
    %         text2 = sprintf('Wcorr p =< %.3f',10^p_val_threshold(index));
    %
    %         lgd = legend([h(1) h(2)],{text1,text2},'Location','northwest');
    %     end
    %
    %     sgtitle('Spearman p =< 0.05 vs Weighted Correlation at equivalent number of replay events detected')
    %     filename = 'Spearman 0.05 comparisions.pdf';
    %     cd ground_truth_original\Figure
    %     saveas(gcf,filename)
    %     cd ..
    %     cd ..
    %     clf
end

end

