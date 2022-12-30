function [] = plot_ground_truth_compare_jump_distance_backup_15_06(folders,option,method)

total_number(1:3) = 0;
for nfolder = 1:length(folders)
    
    cd(folders{nfolder})
    load decoded_replay_events
    
    [~,time_range]=sort_replay_events([],'wcorr');
    
    for event = 1:length(decoded_replay_events(1).replay_events)
        
        event_time = decoded_replay_events(1).replay_events(event).timebins_edges(1);
        
        if event_time >= time_range.pre(1) & event_time <= time_range.pre(2) %If PRE
            
            total_number(1) = total_number(1) + 1;
            
        elseif event_time >= time_range.post(1) & event_time <= time_range.post(2) %If POST
            
            total_number(3) = total_number(3) + 1;
            
        elseif event_time >= time_range.track(1).behaviour(1) & event_time <= time_range.track(1).behaviour(2)
            total_number(2) = total_number(2) + 1;
        elseif event_time >= time_range.track(2).behaviour(1) & event_time <= time_range.track(2).behaviour(2)
            total_number(2) = total_number(2) + 1;
            
            
        end
    end
    cd ..
end

index = [];
epoch_index = [];

for nmethod = 1:length(method)

    cd ground_truth_original
    if strcmp(method{nmethod},'wcorr')
        load log_odd_wcorr_place_bin_circular_shift

    elseif strcmp(method{nmethod},'spearman')
        load log_odd_spearman_place_bin_circular_shift

    elseif strcmp(method{nmethod},'linear')
        load log_odd_linear_place_bin_circular_shift

    elseif strcmp(method{nmethod},'path')
        load log_odd_path_place_bin_circular_shift

    end
    cd ..


    % Get index for PRE,  RUN, POST
    % states = [-1 0 1 2 3 4 5]; % sleepPRE, awakePRE, RUN T1, RUN T2, sleepPOST and awakePOST and multi-track
    % states = [-1 0 1 2 3 4]; % sleepPRE, awakePRE, RUN T1, RUN T2, sleepPOST and awakePOST
    states = [-1 0 1 2];
    jump_threshold = [1 0.6 0.4 0.2];

    for nshuffle = 1:length(jump_threshold)
        for k=1:length(states) % sort track id and behavioral states (pooled from 5 sessions)
            % Find intersect of behavioural state and ripple peak threshold
            state_index = intersect(find(log_odd.behavioural_state==states(k)),...
                find(log_odd.max_jump <= 20*jump_threshold(nshuffle)));

            if k == 1 % PRE
                epoch_index{nmethod}{nshuffle}{1}{1} = intersect(state_index,find(log_odd.track==1));
                epoch_index{nmethod}{nshuffle}{1}{2} = intersect(state_index,find(log_odd.track==2));
            elseif k == 2 % POST
                epoch_index{nmethod}{nshuffle}{3}{1} = intersect(state_index,find(log_odd.track==1));
                epoch_index{nmethod}{nshuffle}{3}{2} = intersect(state_index,find(log_odd.track==2));
            elseif k == 3 % RUN Track 1
                epoch_index{nmethod}{nshuffle}{2}{1} = intersect(state_index,find(log_odd.track==1));
                epoch_index{nmethod}{nshuffle}{2}{2} = intersect(state_index,find(log_odd.track==2));
            elseif k == 4 % RUN Track 2
                epoch_index{nmethod}{nshuffle}{2}{1} = [epoch_index{nmethod}{nshuffle}{2}{1} intersect(state_index,find(log_odd.track==1));];
                epoch_index{nmethod}{nshuffle}{2}{2} = [epoch_index{nmethod}{nshuffle}{2}{2} intersect(state_index,find(log_odd.track==2));];
            end
        end


    end

    % index = epoch_index;

    if strcmp(option,'common')
        for nshuffle = 1:length(jump_threshold)
            data{nmethod}{nshuffle} = log_odd.common_zscored.original(2,:);
            track_id{nmethod}{nshuffle}{1} = find(log_odd.track == 1);
            track_id{nmethod}{nshuffle}{2} = find(log_odd.track == 2);
            event_index{nmethod}{nshuffle} = log_odd.index;
            log_pval{nmethod}{nshuffle} = log_odd.pvalue;
            segment_id{nmethod}{nshuffle} = log_odd.segment_id;

            for session = 1:max(log_odd.experiment)
                session_index{nmethod}{nshuffle}{session} = find(log_odd.experiment == session);
            end
        end

    elseif strcmp(option,'original')
        for nshuffle = 1:length(jump_threshold)
            data{nmethod}{nshuffle} = log_odd.normal_zscored.original(2,:);
            track_id{nmethod}{nshuffle}{1} = find(log_odd.track == 1);
            track_id{nmethod}{nshuffle}{2} = find(log_odd.track == 2);
            event_index{nmethod}{nshuffle} = log_odd.index;
            log_pval{nmethod}{nshuffle} = log_odd.pvalue;
            segment_id{nmethod}{nshuffle} = log_odd_compare{nshuffle}.segment_id;

            for session = 1:max(log_odd.experiment)
                session_index{nmethod}{nshuffle}{session} = find(log_odd.experiment == session);
            end
        end
    end

    %     colour_line2 = {'k--','b--','r--','g--'};
    colour_symbol={'bo','ro','go','ko'};
    Behavioural_epoches = {'PRE','RUN','POST'};
end

p_val_threshold = round(10.^[-3:0.01:log10(1)],4);
p_val_threshold = flip(unique(p_val_threshold));
p_val_threshold(1) = 1;

log_odd_difference = [];
log_odd_difference_CI = [];
percent_sig_events = [];
percent_multi_events = [];

for nmethod = 1:length(method)

    cd ground_truth_original
    if exist('log_odd_difference_jump_distance.mat', 'file') ~= 2
%         load log_odd_difference_single_shuffle
%         % copy place bin shuffle data
%         log_odd_difference1{nmethod}{1} = log_odd_difference{nmethod}{3};
%         log_odd_difference_CI1{nmethod}{1} = log_odd_difference_CI{nmethod}{3};
%         percent_sig_events1{nmethod}{1} = percent_sig_events{nmethod}{3};
%         percent_sig_events_CI1{nmethod}{1} = percent_sig_events_CI{nmethod}{3};
%         percent_multi_events1{nmethod}{1} = percent_multi_events{nmethod}{3};
%         percent_multi_events_CI1{nmethod}{1} = percent_multi_events_CI{nmethod}{3};

        for epoch = 1:3
            for nshuffle = 1:length(log_pval{nmethod})
                %                 multi_event_index = find(log_odd_compare{nshuffle}.multi_track == 1);
                tic
                out = calculate_replay_detection_performance(data{nmethod}{nshuffle},epoch_index{nmethod}{nshuffle}{epoch},segment_id{nmethod}{nshuffle},...
                    log_pval{nmethod}{nshuffle},p_val_threshold,session_index{nmethod}{nshuffle},event_index{nmethod}{nshuffle},total_number(epoch),'bootstrap');
                toc

                log_odd_difference1{nmethod}{nshuffle}{epoch} = out.log_odd_difference;
                log_odd_difference_CI1{nmethod}{nshuffle}{epoch} = prctile(out.log_odd_difference,[2.5 97.5],2);
                percent_sig_events1{nmethod}{nshuffle}{epoch} = out.percent_sig_events;
                percent_sig_events_CI1{nmethod}{nshuffle}{epoch} = prctile(out.percent_sig_events,[2.5 97.5],2);
                percent_multi_events1{nmethod}{nshuffle}{epoch} = out.percent_multi_events;
                percent_multi_events_CI1{nmethod}{nshuffle}{epoch} = prctile(out.percent_multi_events,[2.5 97.5],2);

            end

        end

        log_odd_difference = log_odd_difference1;
        log_odd_difference_CI = log_odd_difference_CI1;
        percent_sig_events = percent_sig_events1;
        percent_sig_events_CI = percent_sig_events_CI1;
        percent_multi_events = percent_multi_events1;
        percent_multi_events_CI = percent_multi_events_CI1;

        save log_odd_difference_jump_distance log_odd_difference log_odd_difference_CI percent_sig_events percent_multi_events percent_sig_events_CI percent_multi_events_CI
        clear log_odd_difference1 log_odd_difference_CI1 percent_sig_events1 percent_sig_events_CI1 percent_multi_events1 percent_multi_events_CI1
        cd ..
    else
        load log_odd_difference_jump_distance
        cd ..
    end


    cd ground_truth_original
    if exist('log_odd_difference_jump_distance_original.mat', 'file') ~= 2
        for epoch = 1:3
            for nshuffle = 1:length(log_pval{nmethod})
                %                 multi_event_index = find(log_odd_compare{nshuffle}.multi_track == 1);
                tic
                out = calculate_replay_detection_performance(data{nmethod}{nshuffle},epoch_index{nmethod}{nshuffle}{epoch},segment_id{nmethod}{nshuffle},...
                    log_pval{nmethod}{nshuffle},p_val_threshold,session_index{nmethod}{nshuffle},event_index{nmethod}{nshuffle},total_number(epoch),'original');

                log_odd_difference_original{nmethod}{nshuffle}{epoch} = out.log_odd_difference;
                percent_sig_events_original{nmethod}{nshuffle}{epoch} = out.percent_sig_events;
                percent_multi_events_original{nmethod}{nshuffle}{epoch} = out.percent_multi_events;
                toc

            end
        end

        save log_odd_difference_jump_distance_original log_odd_difference_original percent_sig_events_original percent_multi_events_original
        cd ..

    else
        load log_odd_difference_jump_distance_original
        cd ..
    end



    p_val_threshold = round(10.^[-3:0.01:log10(1)],4);
    p_val_threshold = flip(unique(p_val_threshold));
    p_val_threshold(1) = 1;

    colour_line= {[127,205,187]/255,[65,182,196]/255,[34,94,168]/255,[37,52,148]/255};
    % colour_line= {[50,136,189]/255,[26,152,80]/255,[244,109,67]/255,[213,62,79]/255};
    p_value_to_plot = round(10.^[-3:0.1:log10(1)],4);
    % p_val_threshold = 0.001:0.0001:0.2;
    p_value_to_plot = flip(unique(p_value_to_plot));
    p_value_to_plot(1) = 1;
    alpha_level = linspace(0.1,0.7,length(p_value_to_plot));
    low_threshold = find(ismembertol(p_val_threshold,p_value_to_plot,eps) == 1);
    jump_distance = {'No Jump distance threshold','Jump distance threshold 0.6','Jump distance threshold 0.4','Jump distance threshold 0.2'};

    fig = figure(1)
    fig.Position = [834 116 850 700];
    for epoch = 1:3
        subplot(2,2,epoch)
        for nshuffle = 1:4
            x = log_odd_difference_original{nmethod}{nshuffle}{epoch}';
            y = percent_sig_events_original{nmethod}{nshuffle}{epoch}';

            for threshold = 1:length(p_value_to_plot)
                if p_value_to_plot(threshold) > 0.0501
                    sc(nshuffle) = scatter(x(low_threshold(threshold)),y(low_threshold(threshold)),'filled','MarkerFaceColor',colour_line{nshuffle},'MarkerFaceAlpha',alpha_level(threshold))
                    hold on
                else
                    sc(nshuffle) = scatter(x(low_threshold(threshold)),y(low_threshold(threshold)),'filled','MarkerFaceColor',colour_line{nshuffle},'MarkerFaceAlpha',alpha_level(threshold),'MarkerEdgeColor',colour_line{nshuffle})
                    hold on
                end
            end
            %             plot(x,y,colour_line{nshuffle});


            xlabel('mean log odd difference')
            ylabel('Proportion of all replay events')
        end

        ylim([0 1])
        title(Behavioural_epoches{epoch});
    end
    legend([sc(1),sc(2),sc(3),sc(4)], {jump_distance{1},jump_distance{2},jump_distance{3},jump_distance{4}})

    cd ground_truth_original\Figure
    filename = sprintf('%s jump distance log odd comparision original.pdf',method{nmethod})
    saveas(gcf,filename)
    clf
    cd ..
    cd ..



    fig = figure(2)
    fig.Position = [834 116 850 700];

    Behavioural_epoches = {'PRE','RUN','POST'};
    low_threshold = find(ismembertol(p_val_threshold,[0.0501 0.02 0.01 0.005 0.002 0.001],eps) == 1);
    colour_line= {[127,205,187]/255,[65,182,196]/255,[34,94,168]/255,[37,52,148]/255};
    alpha_level = linspace(0.1,0.7,length(low_threshold));

    jump_distance = {'No jump distance threshold','Jump distance threshold 0.6','Jump distance threshold 0.4','Jump distance threshold 0.2'};
%     colour_line = {[37,52,148]/255,[34,94,168]/255,[29,145,192]/255,[65,182,196]/255};
    %     colour_line = flip({[12,44,132]/255,'',[29,145,192]/255,'','',[127,205,187]/255});
    
    for epoch = 1:3
        subplot(2,2,epoch)
        for nshuffle = 1:4
            x = mean(log_odd_difference{nmethod}{nshuffle}{epoch},2)';
            y = mean(percent_sig_events{nmethod}{nshuffle}{epoch},2)';

            for threshold = 1:length(low_threshold)
                scatter(x(low_threshold(threshold)),y(low_threshold(threshold)),'filled','MarkerFaceColor',colour_line{nshuffle},'MarkerFaceAlpha',alpha_level(threshold))
                hold on
            end

            UCI = log_odd_difference_CI{nmethod}{nshuffle}{epoch}(:,2)';
            UCI(isnan(x)) = [];
            LCI = log_odd_difference_CI{nmethod}{nshuffle}{epoch}(:,1)';
            LCI(isnan(x)) = [];

            y(isnan(x)) = [];
            x(isnan(x)) = [];
            p(nshuffle) = patch([LCI fliplr(UCI)], [y fliplr(y)], colour_line{nshuffle},'FaceAlpha','0.3','LineStyle','none');
            hold on
            %             plot(x,y,colour_line{nshuffle});

            xlabel('mean log odd difference')
            ylabel('Proportion of all replay events')
        end

        ylim([0 1])
        title(Behavioural_epoches{epoch});
    end

    %     legend([p(1) p(3) p(6)], {'Shuffle 3','Shuffle 2+3','Shuffle 1+2+3'})
    legend([p(1) p(2) p(3) p(4)], {jump_distance{1},jump_distance{2},jump_distance{3},jump_distance{4}},'Position',[0.7 0.3 0.05 0.05])
    %     legend([p(1) p(3) p(6)], {'PRE Place','PRE Place + POST Place','PRE Place + PRE Spike + POST Place'})
    %     legend([p(1),p(2),p(3),p(4),p(5),p(6)], {'Shuffle 3','Shuffle 3+4','Shuffle 2+3','Shuffle 1+3','Shuffle 1+2','Shuffle 1+2+3'})
    cd ground_truth_original\Figure
    filename = sprintf('%s jump distance comparisions CI.pdf',method{nmethod})
    saveas(gcf,filename)
    cd ..
    cd ..
    clf
    %

    
    colour_line= {[127,205,187]/255,[65,182,196]/255,[34,94,168]/255,[37,52,148]/255};
    Behavioural_epoches = {'PRE','RUN','POST'};
    marker_shape = {'o','^','s'};
    
    fig = figure(2)
    fig.Position = [834 116 850 885];
    for nshuffle = 1:4
        for epoch = 1:3
            subplot(2,2,epoch)
            y = mean(percent_multi_events{nmethod}{nshuffle}{epoch}(low_threshold,:),2);
            s(nshuffle) = scatter(p_val_threshold(low_threshold),y,10,marker_shape{epoch},...
                'MarkerFaceColor',colour_line{nshuffle},'MarkerEdgeColor',colour_line{nshuffle},'MarkerFaceAlpha','0.5')
            hold on
            m{nshuffle} = errorbar(p_val_threshold(low_threshold),y,y-percent_multi_events_CI{nmethod}{nshuffle}{epoch}(low_threshold,1),...
                percent_multi_events_CI{nmethod}{nshuffle}{epoch}(low_threshold,2)-y,'.','MarkerFaceColor',colour_line{nshuffle},'MarkerEdgeColor',colour_line{nshuffle})
            m{nshuffle}.Color = colour_line{nshuffle};
            %             hold on
            %             r(nshuffle) = rectangle('Position',pos)

            xlim([0 0.06])
            ylim([0 0.2])
            ylabel('Proportion of multitrack events')
            hold on
            xlabel('Sequenceness p value')
            title(Behavioural_epoches{epoch})
            %             title('Multitrack event proportion vs p value')
        end


        % Sig proportion and log odd at multitrack-corrected p value
        subplot(2,2,4)
        [c index(nshuffle)] = min(abs(mean(percent_multi_events{nmethod}{nshuffle}{epoch},2) - 0.05));

        %         tempt = find(mean(percent_multi_events{nmethod}{nshuffle}{epoch},2) <= 0.05);

        for epoch = 1:3
            x = mean(mean(log_odd_difference{nmethod}{nshuffle}{epoch}(index(nshuffle),:),2),2);
            y = mean(percent_sig_events{nmethod}{nshuffle}{epoch}(index(nshuffle),:),2);
            s(nshuffle) = scatter(x,y,10,marker_shape{epoch},...
                'MarkerFaceColor',colour_line{nshuffle},'MarkerEdgeColor',colour_line{nshuffle},'MarkerFaceAlpha','0.5')
            %             hold on
            %             m{nshuffle} = errorbar(p_val_threshold(low_threshold),y,y-percent_multi_events_CI{nmethod}{nshuffle}{epoch}(low_threshold,1),...
            %                 percent_multi_events_CI{nmethod}{nshuffle}{epoch}(low_threshold,2)-y,'.','MarkerFaceColor',colour_line{nshuffle},'MarkerEdgeColor',colour_line{nshuffle})
            %             m{nshuffle}.Color = colour_line{nshuffle};
            hold on
            error_x = [log_odd_difference_CI{nmethod}{nshuffle}{epoch}(index(nshuffle),1) log_odd_difference_CI{nmethod}{nshuffle}{epoch}(index(nshuffle),2)...
                log_odd_difference_CI{nmethod}{nshuffle}{epoch}(index(nshuffle),2) log_odd_difference_CI{nmethod}{nshuffle}{epoch}(index(nshuffle),1)]
            error_y = [percent_sig_events_CI{nmethod}{nshuffle}{epoch}(index(nshuffle),1) percent_sig_events_CI{nmethod}{nshuffle}{epoch}(index(nshuffle),1)...
                percent_sig_events_CI{nmethod}{nshuffle}{epoch}(index(nshuffle),2) percent_sig_events_CI{nmethod}{nshuffle}{epoch}(index(nshuffle),2)]

            patch(error_x,error_y,colour_line{nshuffle},'EdgeColor',colour_line{nshuffle},'FaceAlpha','0.2','EdgeAlpha','0.2')
            text{nshuffle} = sprintf('%s p =< %.3f (proportion = %.3f)',jump_distance{nshuffle},p_val_threshold(index(nshuffle)),mean(percent_multi_events{nmethod}{nshuffle}{epoch}(index(nshuffle),:)));

            ylabel('Proportion of significant events')
            hold on
            xlabel('mean log odd difference')
            %             title(Behavioural_epoches{epoch})
            title('Equivalent p value when multitrack event proportion = 0.05')
        end
        %         plot([0 0.1],[0 0.3],'r')

        xlim([-1 3])
        ylim([0 0.5])

    end
%     legend([s(1),s(2),s(3),s(4)], {text{1},text{2},text{3},text{4}})
    legend([s(1) s(2) s(3) s(4)], {text{1},text{2},text{3},text{4}})
    %     legend([s(1),s(2),s(3),s(4)], {'spike train circular shift','place field circular shift','place bin circular shift','time bin permutation'},'Position',[0.35 0.2 0.05 0.05])
    %     sgtitle('Multitrack event proportion vs Sequenceness p value')
    cd ground_truth_original\Figure
    filename = sprintf('%s Jump distance multitrack and significant proportion.pdf',method{nmethod})
    saveas(gcf,filename)
    cd ..
    cd ..
    clf

    
%     colour_line = {[127,205,187]/255,'','',[29,145,192]/255,'',[37,52,148]/255};
    fig = figure(3)
    fig.Position = [834 116 860 860];
    for epoch = 1:3
        subplot(2,2,epoch)
        tempt = [mean(percent_sig_events{nmethod}{1}{epoch}(31,:)) mean(percent_sig_events{nmethod}{2}{epoch}(31,:)) mean(percent_sig_events{nmethod}{3}{epoch}(31,:)) ...
            mean(percent_sig_events{nmethod}{4}{epoch}(31,:))];

        [minimum mindex] = min(tempt);

        if epoch == 1
            nedges = -0.5:0.05:0.5;
        elseif epoch == 2
            nedges = 0:0.05:0.8
        elseif epoch == 3
            nedges = 0:0.05:1.8
        end

        for nshuffle = 1:4
            [c index(nshuffle)] = min(abs(mean(percent_sig_events{nmethod}{nshuffle}{epoch},2)-mean(percent_sig_events{nmethod}{mindex}{epoch}(low_threshold(1),:)))); % Get p value index at roughly equivalent propotion of events
            text{nshuffle} = sprintf('%s p =< %.3f (proportion = %.3f)',jump_distance{nshuffle},p_val_threshold(index(nshuffle)),mean(percent_sig_events{nmethod}{nshuffle}{epoch}(index(nshuffle),:)));
            h(nshuffle) = histogram(log_odd_difference{nmethod}{nshuffle}{epoch}(index(nshuffle),:),10,'FaceColor',colour_line{nshuffle},'FaceAlpha',0.3,'Normalization','probability');
            box off
            hold on
        end

        xlabel('mean log odd difference')
        ylabel('Propotion')
        ylim([0 0.6])

        titlename = sprintf('%s',Behavioural_epoches{epoch});
        title(titlename);

        lgd = legend([h(1) h(2) h(3) h(4)],{text{1},text{2},text{3},text{4}},'Location','northwest');
        %        lgd = legend([h(1) h(2) h(3) h(4) h(5) h(6)],{text{1},text{2},text{3},text{4},text{5},text{6}},'Location','northwest');
    end
    titlename = sprintf('%s Jump distance comparision at equivalent proportion',method{nmethod});
    sgtitle(titlename)
    cd ground_truth_original\Figure
    filename = sprintf('%s Jump distance comparision histogram.pdf',method{nmethod});
    saveas(gcf,filename)
    clf
    cd ..
    cd ..

    %% Wcorr ripple comparision
%     shuffle_type = {'POST Place','POST place + POST time','PRE place + POST time','PRE place + POST place','PRE place + PRE spike','PRE place + PRE spike + POST Place'}';
    ripple_type = {'Ripple Threshold 1','Ripple Threshold 3','Ripple Threshold 5','Ripple Threshold 10'}';
%         shuffle_type = {'place field circular shift','spike train circular shift','place bin circular shift','time bin circular shift'}';
%     shuffle_type = {'POST Place','','','PRE place + POST place','','PRE place + PRE spike + POST Place'}';
    p_val_threshold = round(10.^[-3:0.01:-1],4);
    p_val_threshold = flip(unique(p_val_threshold));
    low_threshold = find(ismembertol(p_val_threshold,[0.0501 0.02 0.01 0.005 0.002 0.001],eps) == 1);
    
    ripple_name(1) = {'Ripple threshold'};
    behave_state(1) = {'Behavioural State'};

    Mean_5(1) = {'Mean Log Odd difference at p =< 0.05'};
    Mean_e(1) = {'Mean Log Odd difference at equivalent significant proportion'};
    Mean_m(1) = {'Mean Log Odd difference at equivalent multi-track proportion'};
    
    UCI_5(1) = {'Upper CI at p =< 0.05'};
    UCI_e(1) = {'Upper CI at equivalent significant proportion'};
    UCI_m(1) = {'Upper CI at equivalent multi-track proportion'};

    LCI_5(1) = {'Lower CI at p =< 0.05'};
    LCI_e(1) = {'Lower CI at equivalent significant proportion'};
    LCI_m(1) = {'Lower CI at equivalent multi-track proportion'};

    proportion_5(1) = {'Proportion at p =< 0.05'};
    proportion_e(1) = {'Equivalent propotion'};
    proportion_m(1) = {'Proportion of multi-track event'};

    LCI_proportion_5(1) = {'Lower CI Proportion at p =< 0.05'};
    UCI_proportion_5(1) = {'Upper CI Proportion at p =< 0.05'};    
    LCI_proportion_e(1) = {'Lower CI equivalent significant proportion'};
    UCI_proportion_e(1) = {'Upper CI equivalent significant proportion'};
    LCI_proportion_m(1) = {'Lower CI equivalent multi-track proportion'};
    UCI_proportion_m(1) = {'Upper CI equivalent multi-track proportion'};

    pvalue_e(1) = {'p value at equivalent significant proportion'};
    pvalue_m(1) = {'p value at equivalent multi-track proportion'};

    c = 2;

    for epoch = 1:3
        tempt = [mean(percent_sig_events{nmethod}{2}{epoch}(31,:)) mean(percent_sig_events{nmethod}{3}{epoch}(31,:))...
            mean(percent_sig_events{nmethod}{4}{epoch}(31,:))];
        [minimum mindex] = min(tempt);

        for nshuffle = 1:4
            ripple_name(c) = {ripple_type{nshuffle}};
            Mean_5(c) = {mean(log_odd_difference{nmethod}{nshuffle}{epoch}(31,:))};
            a = prctile(log_odd_difference{nmethod}{nshuffle}{epoch}(31,:),[2.5 97.5]);
            LCI_5(c) = {a(1)};
            UCI_5(c) = {a(2)};
            proportion_5(c) = {mean(percent_sig_events{nmethod}{nshuffle}{epoch}(31,:))};
            b = prctile(percent_sig_events{nmethod}{nshuffle}{epoch}(31,:),[2.5 97.5]);
            LCI_proportion_5(c) = {b(1)};
            UCI_proportion_5(c) = {b(2)};


            [xx index] = min(abs(mean(percent_sig_events{nmethod}{nshuffle}{epoch},2) - mean(percent_sig_events{nmethod}{mindex}{epoch}(31,:)))); % Get p value index at roughly equivalent propotion of events
            Mean_e(c) = {mean(log_odd_difference{nmethod}{nshuffle}{epoch}(index,:))};
            a = prctile(log_odd_difference{nmethod}{nshuffle}{epoch}(index,:),[2.5 97.5]);
            LCI_e(c) = {a(1)};
            UCI_e(c) = {a(1)};
            proportion_e(c) = {mean(percent_sig_events{nmethod}{nshuffle}{epoch}(index,:))};
            b = prctile(percent_sig_events{nmethod}{nshuffle}{epoch}(index,:),[2.5 97.5]);
            LCI_proportion_e(c) = {b(1)};
            UCI_proportion_e(c) = {b(2)};
            pvalue_e(c) = {p_val_threshold(index)};

            [xx index] = min(abs(mean(percent_multi_events{nmethod}{nshuffle}{epoch},2) - 0.05));
            Mean_m(c) = {mean(log_odd_difference{nmethod}{nshuffle}{epoch}(index,:))};
            a = prctile(log_odd_difference{nmethod}{nshuffle}{epoch}(index,:),[2.5 97.5]);
            LCI_m(c) = {a(1)};
            UCI_m(c) = {a(1)};
            proportion_m(c) = {mean(percent_sig_events{nmethod}{nshuffle}{epoch}(index,:))};
            b = prctile(percent_sig_events{nmethod}{nshuffle}{epoch}(index,:),[2.5 97.5]);
            LCI_proportion_m(c) = {b(1)};
            UCI_proportion_m(c) = {b(2)};
            pvalue_m(c) = {p_val_threshold(index)};

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
    cd ground_truth_original\Figure\
    Table = table(behave_state',ripple_name',Mean_5',UCI_5',LCI_5',proportion_5',UCI_proportion_5',LCI_proportion_5',...
        Mean_e',UCI_e',LCI_e',proportion_5',UCI_proportion_5',LCI_proportion_5',pvalue_e',...
        Mean_m',UCI_m',LCI_m',proportion_m',UCI_proportion_m',LCI_proportion_m',pvalue_m');
    
    writetable(Table,'ripple effect table.xlsx')
    cd .. 
    cd ..




end

end
