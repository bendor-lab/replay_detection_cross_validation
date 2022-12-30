function [] = plot_ground_truth_resampled(log_odd,option,method)

% Get index for PRE,  RUN, POST
% states = [-1 0 1 2 3 4 5]; % sleepPRE, awakePRE, RUN T1, RUN T2, sleepPOST and awakePOST and multi-track
% states = [-1 0 1 2 3 4]; % sleepPRE, awakePRE, RUN T1, RUN T2, sleepPOST and awakePOST
states = [-1 0 1 2];
 
epoch_index = [];
% for k=1:length(states)
%     state_index = find(log_odd.behavioural_state==states(k));
%     
%     if states(k) == 1 % On track 1
%         epoch_index{2}{1} = [intersect(state_index,find(log_odd.track==1))]; % Local replay
%     elseif states(k) == 2 % On track 2
%         epoch_index{2}{2} =  [intersect(state_index,find(log_odd.track==2))];
%     elseif states(k) == -1
%         epoch_index{1}{1} = [intersect(state_index,find(log_odd.track==1))];
%         epoch_index{1}{2} = [intersect(state_index,find(log_odd.track==2))];
%     elseif states(k) == 0
%         epoch_index{1}{1} = [epoch_index{1}{1} intersect(state_index,find(log_odd.track==1))];
%         epoch_index{1}{2} = [epoch_index{1}{2} intersect(state_index,find(log_odd.track==2))];
%     elseif states(k) == 3
%         epoch_index{3}{1} = [intersect(state_index,find(log_odd.track==1))];
%         epoch_index{3}{2} = [intersect(state_index,find(log_odd.track==2))];
%     elseif states(k) == 4
%         epoch_index{3}{1} = [epoch_index{3}{1} intersect(state_index,find(log_odd.track==1))];
%         epoch_index{3}{2} = [epoch_index{3}{2} intersect(state_index,find(log_odd.track==2))];
%         %     elseif states(k) == 5
%         %         epoch_index{4} = state_index;
%     end
% end

epoch_index = [];
for k=1:length(states) % sort track id and behavioral states (pooled from 5 sessions)
    state_index = find(log_odd.behavioural_state==states(k));
    
    if k == 1 % PRE
        epoch_index{1}{1} = intersect(state_index,find(log_odd.track==1));
        epoch_index{1}{2} = intersect(state_index,find(log_odd.track==2));
    elseif k == 2 % POST
        epoch_index{3}{1} = intersect(state_index,find(log_odd.track==1));
        epoch_index{3}{2} = intersect(state_index,find(log_odd.track==2));
    elseif k == 3 % RUN Track 1
        epoch_index{2}{1} = intersect(state_index,find(log_odd.track==1)); % Track 1 replay on Track 1 considered RUN Replay            
    elseif k == 4 % RUN Track 2
        epoch_index{2}{2} = intersect(state_index,find(log_odd.track==2));        
    end
end

index = [];

for epoch = 1:3
    tempt1 = datasample(epoch_index{epoch}{1},1000);
    tempt2 = datasample(epoch_index{epoch}{2},1000);
%     index{epoch} = [epoch_index{epoch}{1} epoch_index{epoch}{2}]
    index{epoch} = [ tempt1 tempt2 ]
end



% index = epoch_index;

if strcmp(option,'original')
    data = log_odd.normal_zscored.original(1,:);
elseif strcmp(option,'global remapped')
    data = log_odd.normal_zscored.global_remapped_original(1,:);
end
% data = log_odd.normal_zscored.original(2,:);
% data = log(log_odd.normal.probability_ratio);
% T1_T2_data = log_odd.normal_zscored.original(2,:);


log_pval = log_odd.pvalue;
colour_line= {'k','b','r','g'};
%     colour_line2 = {'k--','b--','r--','g--'};
colour_symbol={'ko','bo','ro','go'};
Behavioural_epoches = {'PRE','RUN','POST','Multi'};

fig = figure(1)
fig.Position = [834 116 850 700];
    
text1 = sprintf('%s',method)
sgtitle(text1);
for epoch = 1:3
    X = [ones(size(index{epoch},2),1) log_pval(index{epoch})'];
    [b,bint,r,rint,stats] = regress(data(index{epoch})',X);
    % the R-square statistic, the F statistic and p value for the full model, and an estimate of the error variance.
    ground_truth(epoch).R_sqaure = stats(1);
    ground_truth(epoch).F_stat = stats(2);
    ground_truth(epoch).p_val = stats(3);
    ground_truth(epoch).error = stats(4);
       
    subplot(2,2,epoch)
%     [histcount_matrix,Xedges,Yedges] = histcounts2(log_pval(index{epoch}),data(index{epoch}),[-3.2:0.1:-1],[-4:0.5:8]);
%     imagesc(histcount_matrix);
%     xticklabels([[Xedges]])
%     xlabel([Xedges])
%     ylabel([Yedges])

    scatter(log_pval(index{epoch}),data(index{epoch}),2,colour_line{epoch},'MarkerEdgeAlpha',0.2)
    hold on
    p = plot(log_pval(index{epoch}),X*b,colour_line{epoch},'LineWidth',2)
    plot([0 -4],[0 0])
    xlim([-3.1 -1])
    xlabel('log p-value')
    ylabel('log odd')
    ylim([-4 8])
    text1 = sprintf('%s - P value %.2d, F Stat %.2d',Behavioural_epoches{epoch},stats(3),stats(2));
    lgd = legend(p,{text1},'Location','southeast');
%     lgd.FontSize = 12;
    legend boxoff
end


% no of significant events vs different mean log odd at different p value threshold

p_val_threshold = flip([-3:0.1:-1.3]);
alpha_level = linspace(0.01,1,length(p_val_threshold));

subplot(2,2,4)
for epoch = 1:3
%     subplot(2,2,epoch)
    for threshold = 1:length(p_val_threshold)
        % find the event index with p value lower than the threshold
%         current_index = intersect(find(log_pval(index{epoch})< p_val_threshold(threshold)),find(log_pval(index{epoch})> p_val_threshold(threshold+1)));
        current_index = find(log_pval(index{epoch})< p_val_threshold(threshold));
        mean_log_odd(threshold) = mean(data(index{epoch}(current_index)));
        mean_log_odd(isnan(mean_log_odd)) = 0;
        SE_log_odd(threshold) = std(data(index{epoch}(current_index)))/sqrt(length((current_index)));
        SE_log_odd(isnan(SE_log_odd)) = 0;
        no_sig_events(threshold) = length(current_index);
        sc = scatter(mean_log_odd(threshold),no_sig_events(threshold),colour_line{epoch},'filled','MarkerFaceAlpha',alpha_level(threshold))
%         plot(mean_log_odd,no_sig_events,colour_line{epoch})
        hold on
    end
    xlabel('log odd')
    ylabel('Number of significant events')
%     lgd.FontSize = 12;
end

cd ground_truth_original

filename = sprintf('%s (%s, resampled) ground truth plot.pdf',method,option)
saveas(gcf,filename)
filename = sprintf('%s (%s, resampled) ground truth plot.pdf',method,option)
saveas(gcf,filename)
cd ..

clf(fig)

end