function plot_ground_truth_compare_test(folders,option,method)

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
    if strcmp(method{nmethod},'wcorr within')
        load log_odd_wcorr
        log_odd_compare{nmethod} = log_odd;

    elseif strcmp(method{nmethod},'spearman')
        load log_odd_spearman
        log_odd_compare{nmethod} = log_odd;

    elseif strcmp(method{nmethod},'linear')
        load log_odd_linear_PRE_place_POST_time
        %         load log_odd_linear_place_bin_circular_shift
        log_odd_compare{nmethod} = log_odd;

    elseif strcmp(method{nmethod},'wcorr across')
        load log_odd_wcorr_normalised_across_tracks
        log_odd_compare{nmethod} = log_odd;

    end
    cd ..

    % Get index for PRE,  RUN, POST
    % states = [-1 0 1 2 3 4 5]; % sleepPRE, awakePRE, RUN T1, RUN T2, sleepPOST and awakePOST and multi-track
    % states = [-1 0 1 2 3 4]; % sleepPRE, awakePRE, RUN T1, RUN T2, sleepPOST and awakePOST
    states = [-1 0 1 2];


    for k=1:length(states) % sort track id and behavioral states (pooled from 5 sessions)
        % Find intersect of behavioural state and ripple peak threshold
        state_index = find(log_odd_compare{nmethod}.behavioural_state==states(k));

        if k == 1 % PRE
            epoch_index{nmethod}{1} = state_index;
        elseif k == 2 % POST
            epoch_index{nmethod}{3} = state_index;
        elseif k == 3 % RUN Track 1
            epoch_index{nmethod}{2} = state_index;
        elseif k == 4 % RUN Track 2
            epoch_index{nmethod}{2} = [epoch_index{nmethod}{2} state_index];
        end
    end

    if strcmp(option,'common')
        for nmethod = 1:length(log_odd_compare)
            data{nmethod} = log_odd_compare{nmethod}.common_zscored.original;
            %             track_id{nmethod}{nshuffle}{1} = find(log_odd_compare{nshuffle}.track == 1);
            %             track_id{nmethod}{nshuffle}{2} = find(log_odd_compare{nshuffle}.track == 2);
            %             event_index{nmethod}{nshuffle} = log_odd_compare{nshuffle}.index;
            log_pval{nmethod} = log_odd_compare{nmethod}.pvalue;
            segment_id{nmethod} = log_odd_compare{nmethod}.segment_id;

            %             for session = 1:max(log_odd_compare{nshuffle}.experiment)
            %                 session_index{nmethod}{nshuffle}{session} = find(log_odd_compare{nshuffle}.experiment == session);
            %             end
        end

    elseif strcmp(option,'original')
        for nmethod = 1:length(log_odd_compare)
            data{nmethod} = log_odd_compare{nmethod}.normal_zscored.original;
            log_pval{nmethod}= log_odd_compare{nmethod}.pvalue;
            %             track_id{nmethod}{nshuffle}{1} = log_odd_compare{nshuffle}.index(find(log_odd_compare{nshuffle}.track == 1));
            %             track_id{nmethod}{nshuffle}{2} = log_odd_compare{nshuffle}.index(find(log_odd_compare{nshuffle}.track == 2));
            %             event_index{nmethod}{nshuffle} = log_odd_compare{nshuffle}.index;
            segment_id{nmethod} = log_odd_compare{nmethod}.segment_id;

            %             for session = 1:max(log_odd_compare{nshuffle}.experiment)
            %                 session_index{nmethod}{nshuffle}{session} = log_odd_compare{nshuffle}.index(find(log_odd_compare{nshuffle}.experiment == session));
            %             end

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
if exist('log_odd_difference_compare_test.mat', 'file') ~= 2


    for epoch = 1:3
        for nmethod = 1:4
            tic
            parfor threshold = 1:length(p_val_threshold)
                for nboot = 1:1000 % Bootstrapping 1000 times

                    %                         resampled_event = epoch_index{nmethod}{nshuffle}{epoch};
                    resampled_event = datasample(epoch_index{nmethod}{epoch},length(epoch_index{nmethod}{epoch}));
                    %                         jump_threshold(nshuffle)

                    %     for track = 1:2
                    this_segment_event = resampled_event(find(segment_id{nmethod}(1,resampled_event) == 1));
                    track_1_index = this_segment_event(find(log_pval{nmethod}(1,this_segment_event) <= p_val_threshold(threshold)));
                    this_segment_event = resampled_event(find(segment_id{nmethod}(1,resampled_event) > 1));
                    track_1_index = [track_1_index this_segment_event(find(log_pval{nmethod}(1,this_segment_event) <= p_val_threshold(threshold)/2))];


                    this_segment_event = resampled_event(find(segment_id{nmethod}(2,resampled_event) == 1));
                    track_2_index = this_segment_event(find(log_pval{nmethod}(2,this_segment_event) <= p_val_threshold(threshold)));
                    this_segment_event = resampled_event(find(segment_id{nmethod}(2,resampled_event) > 1));
                    track_2_index = [track_2_index this_segment_event(find(log_pval{nmethod}(2,this_segment_event) <= p_val_threshold(threshold)/2))];


                    %     end
                    multi_event_number = sum(ismember(track_1_index,track_2_index));
                    multi_event_percent = multi_event_number/(length(track_1_index)+length(track_2_index)-multi_event_number);

                    boot_log_odd_difference(threshold,nboot) = mean(data{nmethod}(1,track_1_index)) ...
                        - mean(data{nmethod}(2,track_2_index));

                    boot_percent_multi_events(threshold,nboot) = multi_event_percent;
                    % Significant event proportion (minus multitrack event number to avoid double counting)
                    boot_percent_sig_events(threshold,nboot) = (length(track_1_index)+length(track_2_index)-multi_event_number)/total_number(epoch);
                end
            end

            log_odd_difference1{nmethod}{epoch} = boot_log_odd_difference;
            log_odd_difference_CI1{nmethod}{epoch} = prctile(boot_log_odd_difference,[2.5 97.5],2);
            percent_sig_events1{nmethod}{epoch} = boot_percent_sig_events;
            percent_sig_events_CI1{nmethod}{epoch} = prctile(boot_percent_sig_events,[2.5 97.5],2);
            percent_multi_events1{nmethod}{epoch} = boot_percent_multi_events;
            percent_multi_events_CI1{nmethod}{epoch} = prctile(boot_percent_multi_events,[2.5 97.5],2);

            toc
        end
    end

    log_odd_difference = log_odd_difference1;
    log_odd_difference_CI = log_odd_difference_CI1;
    percent_sig_events = percent_sig_events1;
    percent_sig_events_CI = percent_sig_events_CI1;
    percent_multi_events = percent_multi_events1;
    percent_multi_events_CI = percent_multi_events_CI1;

    save log_odd_difference_compare_test log_odd_difference log_odd_difference_CI percent_sig_events percent_multi_events percent_sig_events_CI percent_multi_events_CI
    clear log_odd_difference1 log_odd_difference_CI1 percent_sig_events1 percent_sig_events_CI1 percent_multi_events1 percent_multi_events_CI1
    clear boot_log_odd_difference boot_percent_sig_events boot_percent_multi_events
    cd ..
else
    load log_odd_difference_compare_test
    cd ..
end

cd ground_truth_original
if exist('log_odd_difference_compare_test_original.mat', 'file') ~= 2
    for epoch = 1:3
        for nmethod = 1:length(method)
            tic
            for threshold = 1:length(p_val_threshold)
                resampled_event = epoch_index{nmethod}{epoch};
                %                 resampled_event = datasample(epoch_index{nmethod}{epoch},length(epoch_index{nmethod}{epoch}));
                %                         jump_threshold(nshuffle)

                %     for track = 1:2
                this_segment_event = resampled_event(find(segment_id{nmethod}(1,resampled_event) == 1));
                track_1_index = this_segment_event(find(log_pval{nmethod}(1,this_segment_event) <= p_val_threshold(threshold)));
                this_segment_event = resampled_event(find(segment_id{nmethod}(1,resampled_event) > 1));
                track_1_index = [track_1_index this_segment_event(find(log_pval{nmethod}(1,this_segment_event) <= p_val_threshold(threshold)/2))];


                this_segment_event = resampled_event(find(segment_id{nmethod}(2,resampled_event) == 1));
                track_2_index = this_segment_event(find(log_pval{nmethod}(2,this_segment_event) <= p_val_threshold(threshold)));
                this_segment_event = resampled_event(find(segment_id{nmethod}(2,resampled_event) > 1));
                track_2_index = [track_2_index this_segment_event(find(log_pval{nmethod}(2,this_segment_event) <= p_val_threshold(threshold)/2))];


                %     end
                multi_event_number = sum(ismember(track_1_index,track_2_index));
                multi_event_percent = multi_event_number/(length(track_1_index)+length(track_2_index)-multi_event_number);

                boot_log_odd_difference(threshold) = mean(data{nmethod}(1,track_1_index)) ...
                    - mean(data{nmethod}(2,track_2_index));

                boot_percent_multi_events(threshold) = multi_event_percent;
                % Significant event proportion (minus multitrack event number to avoid double counting)
                boot_percent_sig_events(threshold) = (length(track_1_index)+length(track_2_index)-multi_event_number)/total_number(epoch);
            end

            log_odd_difference_original{nmethod}{epoch} = boot_log_odd_difference;
            percent_sig_events_original{nmethod}{epoch} = boot_percent_sig_events;
            percent_multi_events_original{nmethod}{epoch} = boot_percent_multi_events;

            toc
        end
    end

    save log_odd_difference_compare_test_original log_odd_difference_original percent_sig_events_original percent_multi_events_original
    clear boot_log_odd_difference boot_percent_sig_events boot_percent_multi_events
    cd ..
else
    load log_odd_difference_compare_test_original
    cd ..
end


% colour_line = {[69,117,180]/255,[244,165,130]/255,[215,48,39]/255};
 colour_line= {[8,81,156]/255,[67,147,195]/255,[244,109,67]/255,[215,48,39]/255};
% colour_line= {[50,136,189]/255,[26,152,80]/255,[244,109,67]/255,[213,62,79]/255};
p_value_to_plot = round(10.^[-3:0.1:log10(0.2)],4);
% p_val_threshold = 0.001:0.0001:0.2;
p_value_to_plot = flip(unique(p_value_to_plot));
p_value_to_plot(1) = 0.2;
alpha_level = linspace(0.1,0.70,length(p_value_to_plot));
low_threshold = find(ismembertol(p_val_threshold,p_value_to_plot,eps) == 1);

fig = figure(1)
fig.Position = [834 116 850 700];
for epoch = 1:3
    subplot(2,2,epoch)
    for nmethod = 1:4
        x = log_odd_difference_original{nmethod}{epoch}';
        y = percent_sig_events_original{nmethod}{epoch}';

        for threshold = 1:length(p_value_to_plot)
            if p_value_to_plot(threshold) > 0.0501
                sc(nmethod) = scatter(x(low_threshold(threshold)),y(low_threshold(threshold)),'filled','MarkerFaceColor',colour_line{nmethod},...
                    'MarkerFaceAlpha',alpha_level(threshold))
                hold on
            else
                sc(nmethod) = scatter(x(low_threshold(threshold)),y(low_threshold(threshold)),'filled','MarkerFaceColor',colour_line{nmethod},...
                    'MarkerFaceAlpha',alpha_level(threshold),'MarkerEdgeColor',colour_line{nmethod})
                hold on
            end
        end

        xlabel('mean log odd difference')
        ylabel('Proportion of all replay events')
    end

    ylim([0 1])
    title(Behavioural_epoches{epoch});
end
legend([sc(1) sc(2) sc(3),sc(4)], {method{1},method{2},method{3},method{4}})

cd ground_truth_original\Figure
filename = sprintf('%s test comparision original.pdf',method{nmethod})
saveas(gcf,filename)
clf
cd ..
cd ..


fig = figure(2)
fig.Position = [834 116 850 700];
p_val_threshold = round(10.^[-3:0.01:log10(0.2)],4);
p_val_threshold = flip(unique(p_val_threshold));
p_val_threshold(1) = 0.2;
low_threshold = find(ismembertol(p_val_threshold,[0.0501 0.02 0.01 0.005 0.002 0.001],eps) == 1);
 colour_line= {[8,81,156]/255,[67,147,195]/255,[244,109,67]/255,[215,48,39]/255};
alpha_level = linspace(0.2,0.7,length(low_threshold));
%     colour_line = {[37,52,148]/255,[34,94,168]/255,[29,145,192]/255,[65,182,196]/255};
%     colour_line = flip({[12,44,132]/255,'',[29,145,192]/255,'','',[127,205,187]/255});

for epoch = 1:3
    subplot(2,2,epoch)
    for nmethod = 1:4
        x = mean(log_odd_difference{nmethod}{epoch},2)';
        y = mean(percent_sig_events{nmethod}{epoch},2)';

        for threshold = 1:length(low_threshold)
            scatter(x(low_threshold(threshold)),y(low_threshold(threshold)),'filled','MarkerFaceColor',colour_line{nmethod},...
                'MarkerFaceAlpha',alpha_level(threshold),'MarkerEdgeColor',colour_line{nmethod})
            hold on
        end

        UCI = log_odd_difference_CI{nmethod}{epoch}(:,2)';
        UCI(isnan(x)) = [];
        LCI = log_odd_difference_CI{nmethod}{epoch}(:,1)';
        LCI(isnan(x)) = [];

        y(isnan(x)) = [];
        x(isnan(x)) = [];
        p(nmethod) = patch([LCI fliplr(UCI)], [y fliplr(y)], colour_line{nmethod},'FaceAlpha','0.3','LineStyle','none');
        hold on
        %             plot(x,y,colour_line{nshuffle});

        xlabel('mean log odd difference')
        ylabel('Proportion of all replay events')
    end

    ylim([0 0.7])
    title(Behavioural_epoches{epoch});
end

%     legend([p(1) p(3) p(6)], {'Shuffle 3','Shuffle 2+3','Shuffle 1+2+3'})
% legend([p(1) p(2) p(3)], {method{1},method{2},method{3}},'Position',[0.7 0.3 0.05 0.05])
legend([p(1) p(2) p(3) p(4)], {method{1},method{2},method{3},method{4}})
%     legend([p(1) p(3) p(6)], {'PRE Place','PRE Place + POST Place','PRE Place + PRE Spike + POST Place'})
%     legend([p(1),p(2),p(3),p(4),p(5),p(6)], {'Shuffle 3','Shuffle 3+4','Shuffle 2+3','Shuffle 1+3','Shuffle 1+2','Shuffle 1+2+3'})
cd ground_truth_original\Figure
filename = sprintf('%s test comparisions CI.pdf',method{nmethod})
saveas(gcf,filename)
cd ..
cd ..
clf
%

%     nshuffle
 colour_line= {[8,81,156]/255,[67,147,195]/255,[244,109,67]/255,[215,48,39]/255};
Behavioural_epoches = {'PRE','RUN','POST'};
marker_shape = {'o','^','s'};

fig = figure(2)
fig.Position = [834 116 850 885];
for nmethod = 1:4
    for epoch = 1:3
        subplot(2,2,epoch)
        y = mean(percent_multi_events{nmethod}{epoch}(low_threshold,:),2);
        s(nmethod) = scatter(p_val_threshold(low_threshold),y,10,marker_shape{epoch},...
            'MarkerFaceColor',colour_line{nmethod},'MarkerEdgeColor',colour_line{nmethod},'MarkerFaceAlpha','0.5')
        hold on
        m{nmethod} = errorbar(p_val_threshold(low_threshold),y,y-percent_multi_events_CI{nmethod}{epoch}(low_threshold,1),...
            percent_multi_events_CI{nmethod}{epoch}(low_threshold,2)-y,'.','MarkerFaceColor',colour_line{nmethod},'MarkerEdgeColor',colour_line{nmethod})
        m{nmethod}.Color = colour_line{nmethod};
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

    for epoch = 1:3
        [c index(nmethod)] = min(abs(mean(percent_multi_events{nmethod}{epoch},2) - 0.05));
        x = mean(mean(log_odd_difference{nmethod}{epoch}(index(nmethod),:),2),2);
        y = mean(percent_sig_events{nmethod}{epoch}(index(nmethod),:),2);
        s(nmethod) = scatter(x,y,10,marker_shape{epoch},...
            'MarkerFaceColor',colour_line{nmethod},'MarkerEdgeColor',colour_line{nmethod},'MarkerFaceAlpha','0.5')
        %             hold on
        %             m{nshuffle} = errorbar(p_val_threshold(low_threshold),y,y-percent_multi_events_CI{nmethod}{nshuffle}{epoch}(low_threshold,1),...
        %                 percent_multi_events_CI{nmethod}{nshuffle}{epoch}(low_threshold,2)-y,'.','MarkerFaceColor',colour_line{nshuffle},'MarkerEdgeColor',colour_line{nshuffle})
        %             m{nshuffle}.Color = colour_line{nshuffle};
        hold on
        error_x = [log_odd_difference_CI{nmethod}{epoch}(index(nmethod),1) log_odd_difference_CI{nmethod}{epoch}(index(nmethod),2)...
            log_odd_difference_CI{nmethod}{epoch}(index(nmethod),2) log_odd_difference_CI{nmethod}{epoch}(index(nmethod),1)]
        error_y = [percent_sig_events_CI{nmethod}{epoch}(index(nmethod),1) percent_sig_events_CI{nmethod}{epoch}(index(nmethod),1)...
            percent_sig_events_CI{nmethod}{epoch}(index(nmethod),2) percent_sig_events_CI{nmethod}{epoch}(index(nmethod),2)]

        patch(error_x,error_y,colour_line{nmethod},'EdgeColor',colour_line{nmethod},'FaceAlpha','0.2','EdgeAlpha','0.2')
        text{nmethod} = sprintf('%s p =< %.3f (proportion = %.3f)',method{nmethod},p_val_threshold(index(nmethod)),mean(percent_multi_events{nmethod}{epoch}(index(nmethod),:)));

        ylabel('Proportion of significant events')
        hold on
        xlabel('mean log odd difference')
        %             title(Behavioural_epoches{epoch})
        title('Equivalent p value when multitrack event proportion = 0.05')
    end
    %         plot([0 0.1],[0 0.3],'r')

    xlim([-1 3])
    ylim([0 0.6])
end

%     legend([s(1),s(2),s(3),s(4)], {text{1},text{2},text{3},text{4}})
legend([s(1) s(2) s(3) s(4)], {text{1},text{2},text{3},text{4}})
%     legend([s(1),s(2),s(3),s(4)], {'spike train circular shift','place field circular shift','place bin circular shift','time bin permutation'},'Position',[0.35 0.2 0.05 0.05])
%     sgtitle('Multitrack event proportion vs Sequenceness p value')
cd ground_truth_original\Figure
filename = sprintf('%s test multitrack and significant proportion.pdf',method{nmethod})
saveas(gcf,filename)
cd ..
cd ..
clf


%     colour_line = {[127,205,187]/255,'','',[29,145,192]/255,'',[37,52,148]/255};
fig = figure(3)
fig.Position = [834 116 860 860];
for epoch = 1:3
    subplot(2,2,epoch)
    tempt = [mean(percent_sig_events{1}{epoch}(low_threshold(1),:)) mean(percent_sig_events{2}{epoch}(low_threshold(1),:)) mean(percent_sig_events{3}{epoch}(low_threshold(1),:))];

    [minimum mindex] = min(tempt);
    if epoch == 1
        nedges = -0.5:0.05:0.5;
    elseif epoch == 2
        nedges = 0:0.05:0.8
    elseif epoch == 3
        nedges = 0:0.05:1.8
    end

    for nmethod = 1:3
        [c index(nmethod)] = min(abs(mean(percent_sig_events{nmethod}{epoch},2)-mean(percent_sig_events{mindex}{epoch}(low_threshold(1),:)))); % Get p value index at roughly equivalent propotion of events
        text{nmethod} = sprintf('%s p =< %.3f (proportion = %.3f)',method{nmethod},p_val_threshold(index(nmethod)),mean(percent_sig_events{nmethod}{epoch}(index(nmethod),:)));
        h(nmethod) = histogram(log_odd_difference{nmethod}{epoch}(index(nmethod),:),10,'FaceColor',colour_line{nmethod},'FaceAlpha',0.3,'Normalization','probability');
        box off
        hold on
    end

    xlabel('mean log odd difference')
    ylabel('Propotion')
    ylim([0 0.6])

    titlename = sprintf('%s',Behavioural_epoches{epoch});
    title(titlename);

    lgd = legend([h(1) h(2) h(3)],{text{1},text{2},text{3}},'Location','northwest');
    %        lgd = legend([h(1) h(2) h(3) h(4) h(5) h(6)],{text{1},text{2},text{3},text{4},text{5},text{6}},'Location','northwest');
end
titlename = sprintf('%s 3 methods comparision at equivalent proportion',method{nmethod});
sgtitle(titlename)
cd ground_truth_original\Figure
filename = sprintf('test comparision histogram.pdf',method{nmethod});
saveas(gcf,filename)
clf
cd ..
cd ..


% Compare at equivalent proportion (the method with the lowest proportion of events detected)
method_type = {'Weighted Correlation','Spearman Correlation','Linear Fitting'}';
%     colour_line = {'k','g','c','b','m','r'};
index = [];
% colour_line= {'b','r','g','k'};
% colour_line= {[50,136,189]/255,[26,152,80]/255,[213,62,79]/255};
fig = figure(2)
fig.Position = [834 116 860 860];
for epoch = 1:3
    subplot(2,2,epoch)

    for nmethod = 1:length(method_type)

        text{nmethod} = sprintf('%s p =< %.3f (proportion = %.3f)',method_type{nmethod},p_val_threshold(low_threshold(1)),percent_sig_events{nmethod}{epoch}(1));
        h(nmethod) = histogram(log_odd_difference{nmethod}{epoch}(low_threshold(1),:),10,'FaceColor',colour_line{nmethod},'FaceAlpha',0.3,'Normalization','probability');
        box off
        hold on

    end

    xlabel('mean log odd difference')
    ylabel('Propotion')
    ylim([0 0.6])

    titlename = sprintf('%s',Behavioural_epoches{epoch});
    title(titlename);

    lgd = legend([h(1) h(2) h(3)],{text{1},text{2},text{3}},'Location','northwest');
    %        lgd = legend([h(1) h(2) h(3) h(4) h(5) h(6)],{text{1},text{2},text{3},text{4},text{5},text{6}},'Location','northwest');
end
titlename = 'methods comparision at equivalent p value';
sgtitle(titlename)
cd ground_truth_original\Figure
filename = 'test comparision 0.05 histogram.pdf';
saveas(gcf,filename)
clf
cd ..
cd ..

% %% 3 Methods comparision
% %     shuffle_type = {'POST Place','POST place + POST time','PRE place + POST time','PRE place + POST place','PRE place + PRE spike','PRE place + PRE spike + POST Place'}';
% % ripple_type = {'Ripple Threshold 1','Ripple Threshold 3','Ripple Threshold 5','Ripple Threshold 10'}';
% method_type = {'Weighted correlation','Spearman correlation','Linear fitting'};
% %         shuffle_type = {'place field circular shift','spike train circular shift','place bin circular shift','time bin circular shift'}';
% %     shuffle_type = {'POST Place','','','PRE place + POST place','','PRE place + PRE spike + POST Place'}';
% 
% p_val_threshold = round(10.^[-3:0.01:log10(0.2)],4);
% p_val_threshold = flip(unique(p_val_threshold));
% p_val_threshold(1) = 0.2;
% low_threshold = find(ismembertol(p_val_threshold,[0.0501 0.02 0.01 0.005 0.002 0.001],eps) == 1);
% 
% method_name(1) = {'Methods'};
% behave_state(1) = {'Behavioural State'};
% 
% Mean_5(1) = {'Mean Log Odd difference at p =< 0.05'};
% Mean_e(1) = {'Mean Log Odd difference at equivalent significant proportion'};
% Mean_m(1) = {'Mean Log Odd difference at equivalent multi-track proportion'};
% 
% UCI_5(1) = {'Upper CI at p =< 0.05'};
% UCI_e(1) = {'Upper CI at equivalent significant proportion'};
% UCI_m(1) = {'Upper CI at equivalent multi-track proportion'};
% 
% LCI_5(1) = {'Lower CI at p =< 0.05'};
% LCI_e(1) = {'Lower CI at equivalent significant proportion'};
% LCI_m(1) = {'Lower CI at equivalent multi-track proportion'};
% 
% proportion_5(1) = {'Proportion at p =< 0.05'};
% proportion_e(1) = {'Equivalent propotion'};
% proportion_em(1) = {'Equivalent multi-track propotion'};
% proportion_m(1) = {'Proportion of multi-track event'};
% 
% LCI_proportion_5(1) = {'Lower CI Proportion at p =< 0.05'};
% UCI_proportion_5(1) = {'Upper CI Proportion at p =< 0.05'};
% LCI_proportion_e(1) = {'Lower CI equivalent significant proportion'};
% UCI_proportion_e(1) = {'Upper CI equivalent significant proportion'};
% LCI_proportion_em(1) = {'Lower CI equivalent multi-track proportion'};
% UCI_proportion_em(1) = {'Upper CI equivalent multi-track proportion'};
% LCI_proportion_m(1) = {'Lower CI multi-track proportion'};
% UCI_proportion_m(1) = {'Upper CI multi-track proportion'};
% 
% pvalue_e(1) = {'p value at equivalent significant proportion'};
% pvalue_m(1) = {'p value at equivalent multi-track proportion'};
% 
% c = 2;
% for epoch = 1:3
%     tempt = [mean(percent_sig_events{1}{epoch}(low_threshold(1),:)) mean(percent_sig_events{2}{epoch}(low_threshold(1),:))...
%         mean(percent_sig_events{3}{epoch}(low_threshold(1),:))];
%     [minimum mindex] = min(tempt);
% 
%     for nmethod = 1:3
%         method_name(c) = {method_type{nmethod}};
%         Mean_5(c) = {mean(log_odd_difference{nmethod}{epoch}(low_threshold(1),:))};
%         a = prctile(log_odd_difference{nmethod}{epoch}(low_threshold(1),:),[2.5 97.5]);
%         LCI_5(c) = {a(1)};
%         UCI_5(c) = {a(2)};
%         proportion_5(c) = {mean(percent_sig_events{nmethod}{epoch}(low_threshold(1),:))};
%         b = prctile(percent_sig_events{nmethod}{epoch}(low_threshold(1),:),[2.5 97.5]);
%         LCI_proportion_5(c) = {b(1)};
%         UCI_proportion_5(c) = {b(2)};
% 
% 
%         [xx index] = min(abs(mean(percent_sig_events{nmethod}{epoch},2) - mean(percent_sig_events{mindex}{epoch}(low_threshold(1),:)))); % Get p value index at roughly equivalent propotion of events
%         Mean_e(c) = {mean(log_odd_difference{nmethod}{epoch}(index,:))};
%         a = prctile(log_odd_difference{nmethod}{epoch}(index,:),[2.5 97.5]);
%         LCI_e(c) = {a(1)};
%         UCI_e(c) = {a(1)};
%         proportion_e(c) = {mean(percent_sig_events{nmethod}{epoch}(index,:))};
%         b = prctile(percent_sig_events{nmethod}{epoch}(index,:),[2.5 97.5]);
%         LCI_proportion_e(c) = {b(1)};
%         UCI_proportion_e(c) = {b(2)};
%         pvalue_e(c) = {p_val_threshold(index)};
% 
%         [xx index] = min(abs(mean(percent_multi_events{nmethod}{epoch},2) - 0.05));
%         Mean_m(c) = {mean(log_odd_difference{nmethod}{epoch}(index,:))};
%         a = prctile(log_odd_difference{nmethod}{epoch}(index,:),[2.5 97.5]);
%         LCI_m(c) = {a(1)};
%         UCI_m(c) = {a(1)};
%         
%         proportion_em(c) = {mean(percent_sig_events{nmethod}{epoch}(index,:))};
%         b = prctile(percent_sig_events{nmethod}{epoch}(index,:),[2.5 97.5]);
%         LCI_proportion_em(c) = {b(1)};
%         UCI_proportion_em(c) = {b(2)};
% 
%         proportion_m(c) = {mean(percent_multi_events{nmethod}{epoch}(index,:))};
%         b = prctile(percent_multi_events{nmethod}{epoch}(index,:),[2.5 97.5]);
%         LCI_proportion_m(c) = {b(1)};
%         UCI_proportion_m(c) = {b(2)};
% 
%         pvalue_m(c) = {p_val_threshold(index)};
% 
%         if epoch == 1
%             behave_state(c) = {'PRE'};
%         elseif epoch == 2
% 
%             behave_state(c) = {'RUN'};
%         elseif epoch == 3
% 
%             behave_state(c) = {'POST'};
%         end
%         c = c+1;
%     end
% 
% 
% end
% cd ground_truth_original\Figure\
% Table = table(behave_state',method_name',Mean_5',UCI_5',LCI_5',proportion_5',UCI_proportion_5',LCI_proportion_5',...
%     Mean_e',UCI_e',LCI_e',proportion_5',UCI_proportion_5',LCI_proportion_5',pvalue_e',...
%     Mean_m',UCI_m',LCI_m',proportion_em',UCI_proportion_em',LCI_proportion_em',proportion_m',UCI_proportion_m',LCI_proportion_m',pvalue_m');
% 
% writetable(Table,'3 methods compare table.xlsx')
% cd ..
% cd ..
% 
% 
% % Table for original data
% c = 2;
% for epoch = 1:3
%     tempt = [percent_sig_events_original{1}{epoch}(low_threshold(1)) percent_sig_events_original{2}{epoch}(low_threshold(1))...
%         percent_sig_events_original{3}{epoch}(low_threshold(1))];
%     [minimum mindex] = min(tempt);
% 
%     for nmethod = 1:3
%         method_name(c) = {method_type{nmethod}};
%         Mean_5(c) = {log_odd_difference_original{nmethod}{epoch}(low_threshold(1))};
%         proportion_5(c) = {percent_sig_events_original{nmethod}{epoch}(low_threshold(1))};
% 
% 
%         [xx index] = min(abs(percent_sig_events_original{nmethod}{epoch} - percent_sig_events_original{mindex}{epoch}(low_threshold(1)))); % Get p value index at roughly equivalent propotion of events
%         Mean_e(c) = {log_odd_difference_original{nmethod}{epoch}(index)};
%         proportion_e(c) = {percent_sig_events_original{nmethod}{epoch}(index)};
%         pvalue_e(c) = {p_val_threshold(index)};
% 
%         [xx index] = min(abs(percent_multi_events_original{nmethod}{epoch} - 0.05));
%         Mean_m(c) = {log_odd_difference_original{nmethod}{epoch}(index)};
%         proportion_m(c) = {percent_multi_events_original{nmethod}{epoch}(index)};
%         proportion_em(c) = {percent_sig_events_original{nmethod}{epoch}(index)};
%         pvalue_m(c) = {p_val_threshold(index)};
% 
%         if epoch == 1
%             behave_state(c) = {'PRE'};
%         elseif epoch == 2
% 
%             behave_state(c) = {'RUN'};
%         elseif epoch == 3
% 
%             behave_state(c) = {'POST'};
%         end
%         c = c+1;
%     end
% end
% cd ground_truth_original\Figure\
% Table = table(behave_state',method_name',Mean_5',proportion_5',...
%     Mean_e',proportion_e',pvalue_e',...
%     Mean_m',proportion_em',proportion_m',pvalue_m');
% 
% writetable(Table,'3 methods compare original table.xlsx')
% cd ..
% cd ..

end





