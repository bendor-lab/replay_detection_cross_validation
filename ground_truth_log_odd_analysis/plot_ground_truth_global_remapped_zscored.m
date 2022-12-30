function [] = plot_ground_truth_global_remapped_zscored(folders,method,option)

total_number{1}(1:3) = 0;
total_number{2}(1:3) = 0;

for nfolder = 1:5
    for nshuffle = 1:2
        
        cd(folders{nfolder})
        
        if nshuffle == 2
            cd 100_global_remapped_shuffles
            number_of_global_remapped_shuffles = length(dir('shuffle_*'));
            cd ..
        elseif nshuffle == 1 % Original data, no shuffle
            number_of_global_remapped_shuffles = 1;
        end
        
        load decoded_replay_events
        
        [~,time_range]=sort_replay_events([],'wcorr');
        
        for event = 1:length(decoded_replay_events(1).replay_events)
            
            event_time = decoded_replay_events(1).replay_events(event).timebins_edges(1);
            
            if event_time <= time_range.pre(2) %If PRE
                
                total_number{nshuffle}(1) = total_number{nshuffle}(1) + 1*number_of_global_remapped_shuffles;
                
            elseif event_time >= time_range.post(1) %If POST
                
                total_number{nshuffle}(3) = total_number{nshuffle}(3) + 1*number_of_global_remapped_shuffles;
                
            elseif event_time >= time_range.track(1).behaviour(1) & event_time <= time_range.track(1).behaviour(2)
                total_number{nshuffle}(2) = total_number{nshuffle}(2) + 1*number_of_global_remapped_shuffles;
            elseif event_time >= time_range.track(2).behaviour(1) & event_time <= time_range.track(2).behaviour(2)
                total_number{nshuffle}(2) = total_number{nshuffle}(2) + 1*number_of_global_remapped_shuffles;
                
                
            end
        end
        cd ..
    end
end

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
    index = [];
    states = [-1 0 1 2];

    for nshuffle = 1:2
        epoch_index = [];
        for k=1:length(states) % sort track id and behavioral states (pooled from 5 sessions)
            state_index = find(log_odd_compare{nshuffle}.behavioural_state==states(k));
            
            if k == 1 % PRE
                epoch_index{1}{1} = intersect(state_index,find(log_odd_compare{nshuffle}.track==1));
                epoch_index{1}{2} = intersect(state_index,find(log_odd_compare{nshuffle}.track==2));
            elseif k == 2 % POST
                epoch_index{3}{1} = intersect(state_index,find(log_odd_compare{nshuffle}.track==1));
                epoch_index{3}{2} = intersect(state_index,find(log_odd_compare{nshuffle}.track==2));
            elseif k == 3 % RUN Track 1
                epoch_index{2}{1} = intersect(state_index,find(log_odd_compare{nshuffle}.track==1)); % Track 1 replay on Track 1 considered RUN Replay
            elseif k == 4 % RUN Track 2
                epoch_index{2}{2} = intersect(state_index,find(log_odd_compare{nshuffle}.track==2));
            end
            
        end
        
        for epoch = 1:3
            %         tempt1 = datasample(epoch_index{epoch}{1},1000);
            %         tempt2 = datasample(epoch_index{epoch}{2},1000);
            index{nmethod}{nshuffle}{epoch} = [epoch_index{epoch}{1} epoch_index{epoch}{2}]
            %         index{nmethod}{epoch} = [ tempt1 tempt2 ]
        end
        
    end


    % index = epoch_index;
    if strcmp(option,'original')
        data{nmethod}{1} = log_odd_compare{1}.normal_zscored.original(1,:);
        data{nmethod}{2} = log_odd_compare{2}.normal_zscored.global_remapped_original(1,:);
    elseif strcmp(option,'common')
        data{nmethod}{1} = log_odd_compare{1}.common_zscored.original(1,:);
        data{nmethod}{2} = log_odd_compare{2}.common_zscored.global_remapped_original(1,:);
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
filename = sprintf('%s (%s) Orignal vs shuffle zscored.pdf',method{nmethod},option)
sgtitle(filename);
p_val_threshold = flip([-3:0.1:-1.3]);
alpha_level = linspace(0.01,1,length(p_val_threshold));

for epoch = 1:3
    subplot(2,2,epoch)

    for nshuffle = 1:2

        for threshold = 1:length(p_val_threshold)
            % find the event index with p value lower than the threshold
            %         current_index = intersect(find(log_pval(index{epoch})< p_val_threshold(threshold)),find(log_pval(index{epoch})> p_val_threshold(threshold+1)));
            current_index = find(log_pval{nmethod}{nshuffle}(index{nmethod}{nshuffle}{epoch})< p_val_threshold(threshold));

            % Original
            mean_log_odd{1}{nmethod}{epoch}{nshuffle}(threshold) = mean(data{nmethod}{nshuffle}(index{nmethod}{nshuffle}{epoch}(current_index)));
            mean_log_odd{1}{nmethod}{epoch}{nshuffle}(isnan(mean_log_odd{1}{nmethod}{epoch}{nshuffle})) = 0;
            SE_log_odd{1}{nmethod}{epoch}{nshuffle}(threshold) = std(data{nmethod}{nshuffle}(index{nmethod}{nshuffle}{epoch}(current_index)))/sqrt(length((current_index)));
            SE_log_odd{1}{nmethod}{epoch}{nshuffle}(isnan(SE_log_odd{1}{nmethod}{epoch}{nshuffle})) = 0;
            
            % Zscore normalised relative to the shuffled distribution
            zdata = (data{nmethod}{nshuffle}(index{nmethod}{nshuffle}{epoch}(current_index)) - mean(data{nmethod}{2})) ./ std(data{nmethod}{2});

            mean_log_odd{2}{nmethod}{epoch}{nshuffle}(threshold) = mean(zdata);
            mean_log_odd{2}{nmethod}{epoch}{nshuffle}(isnan(mean_log_odd{2}{nmethod}{epoch}{nshuffle})) = 0;
            SE_log_odd{2}{nmethod}{epoch}{nshuffle}(threshold) = std(zdata)/sqrt(length((current_index)));
            SE_log_odd{2}{nmethod}{epoch}{nshuffle}(isnan(SE_log_odd{2}{nmethod}{epoch}{nshuffle})) = 0;


            %             no_sig_events(threshold) = length(current_index)/length(data{nmethod});
            percent_sig_events{nmethod}{epoch}{nshuffle}(threshold) = length(current_index)/total_number{nshuffle}(epoch);
                %             sc(nshuffle) = scatter(mean_log_odd{nmethod}{nshuffle}(threshold),percent_sig_events{nmethod}{nshuffle}(threshold),colour_line{nshuffle},'filled','MarkerFaceAlpha',alpha_level(threshold))
                %         plot(mean_log_odd,no_sig_events,colour_line{epoch})
                
                
                %             hold on
            end

            
            %     lgd.FontSize = 12;
            
            
            %     legend([sc(1),sc(2),sc(3),sc(4)], {'wcorr','spearman','linear','path'})
            %         title(Behavioural_epoches{epoch});

    end
    
    
    %     ft = fittype( 'poly2' );
    %     opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    %     opts.Display = 'Off';
    %     opts.StartPoint = [1 1];
    %
    %     fit_curve{nmethod} = fit(mean_log_odd{nmethod}',no_sig_events{nmethod}',ft,opts)
    %     plot(fit_curve{nmethod},colour_line{nmethod})
%     hold on
%     xlabel('mean log odd')
%     ylabel('Proportion of all replay events')
    %     lgd.FontSize = 12;
    
%     legend([sc(1),sc(2)], {'Original','Cell ID Shuffled'})
end


for epoch = 1:3
    subplot(2,2,epoch)
    for threshold = 1:length(p_val_threshold)
        for n = 1:4
            
            if n == 1
                sc(n) = scatter(mean_log_odd{1}{nmethod}{epoch}{1}(threshold),percent_sig_events{nmethod}{epoch}{1}(threshold),colour_line{n},'filled','MarkerFaceAlpha',alpha_level(threshold))
                
            elseif n == 2
                sc(n) = scatter(mean_log_odd{2}{nmethod}{epoch}{1}(threshold),percent_sig_events{nmethod}{epoch}{1}(threshold),colour_line{n},'filled','MarkerFaceAlpha',alpha_level(threshold))

            elseif n == 3
                sc(n) = scatter(mean_log_odd{1}{nmethod}{epoch}{2}(threshold),percent_sig_events{nmethod}{epoch}{2}(threshold),colour_line{n},'filled','MarkerFaceAlpha',alpha_level(threshold))
            elseif n == 4
                sc(n) = scatter(mean_log_odd{2}{nmethod}{epoch}{2}(threshold),percent_sig_events{nmethod}{epoch}{2}(threshold),colour_line{n},'filled','MarkerFaceAlpha',alpha_level(threshold))
                           

            end
            
            hold on
            xlabel('mean log odd')
            ylabel('Proportion of all replay events')
        end
    end
    title(Behavioural_epoches{epoch});
    if epoch ==1
        legend([sc(1),sc(2),sc(3),sc(4)], {'Original','Zscored relative to the shuffle','Cell ID shuffle','Cell ID shuffle zscored'})
    end
end

cd ground_truth_original\Figure

filename = sprintf('%s (%s) Orignal vs shuffle zscored.pdf',method{nmethod},option)
saveas(gcf,filename)
% filename = sprintf('%s (%s) ground truth plot.pdf',method,option)
% saveas(gcf,filename)
cd ..
cd ..

clf(fig)

end
end

