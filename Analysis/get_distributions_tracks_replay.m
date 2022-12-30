function get_distributions_tracks_replay(varargin)

available_controls= {'fixed_rate','fixed_spike_count','global_remapped','rate_remapped'};
available_scoring= {'wcorr','spearman'};
available_methods= {'TRACK_PAIRS','ONE_TRACK'};

p= inputParser;
addParameter(p,'scoring','wcorr',@(x) ismember(available_scoring,x));
addParameter(p,'control',[],@(x) ismember(available_controls,x));
addParameter(p,'method','TRACK_PAIRS',@(x) ismember(available_methods,x));
parse(p,varargin{:});

cd('X:\BendorLab\Drobo\Neural and Behavioural Data\Rate remapping\Data')

filename= ['rate_remapping_analysis_' p.Results.method '_' p.Results.scoring];
if isempty(p.Results.control)
    load([filename '.mat']);
else
    switch p.Results.control
        case 'fixed_rate'
            load(['.\CONTROLS\' p.Results.control '\' filename '_FIXED.mat']);
        otherwise
            load(['.\CONTROLS\' p.Results.control '\' filename '.mat']);
    end
end
    num_exp= length(unique(remapping(6).experiment));

    figure('Color','w');
    % mean FR track
    subplot(2*num_exp,2*num_exp+2,  [1,(num_exp+1)+(num_exp-1)*2*(num_exp+1)-1]);%(num_exp+1)+2*(num_exp+1).* 3]); 
    [n_T1,edgesT1]= histcounts(log(remapping_raw(6).mean_rate_T1),'Normalization','probability');
    [n_T2]= histcounts(log(remapping_raw(6).mean_rate_T2),edgesT1,'Normalization','probability');
    ctrs_T1= edgesT1(1:end-1)+ 0.5*mean(diff(edgesT1));
    bar(ctrs_T1,n_T1,'FaceColor',[1 0.2 0.2],'FaceAlpha',0.5,'EdgeAlpha',0); hold on;
    bar(ctrs_T1,n_T2,'FaceColor',[0.2 0.3 1],'FaceAlpha',0.5,'EdgeAlpha',0);
    xlabel('log mean FR (Hz)'); xticks([min(edgesT1)  max(edgesT1)]);
    xticklabels(exp([min(edgesT1)  max(edgesT1)])); xtickformat('%0.1f');
    yticks([]); ylabel('probability');
    set(gca,'Box','off','TickDir','out'); 
    legend({'Track 1','Track2'},'Box','off')
    title('Mean FR on Track')
    [~,ks_stats.MeanFrTrack_all] = kstest2(remapping_raw(6).mean_rate_T1,remapping_raw(6).mean_rate_T2);
    for this_sess=1:num_exp
        subplot(2*num_exp,2*num_exp+2, [(num_exp+1)+ (this_sess-1)*2*(num_exp+1)]);
        n_T1= histcounts(log(remapping_raw(6).mean_rate_T1(remapping_raw(6).experiment == this_sess)),edgesT1,'Normalization','probability');
        n_T2= histcounts(log(remapping_raw(6).mean_rate_T2(remapping_raw(6).experiment == this_sess)),edgesT1,'Normalization','probability');
        bar(ctrs_T1,n_T1,'FaceColor',[1 0.2 0.2],'FaceAlpha',0.5,'EdgeAlpha',0); hold on;
        bar(ctrs_T1,n_T2,'FaceColor',[0.2 0.3 1],'FaceAlpha',0.5,'EdgeAlpha',0);
        set(gca,'Box','off','TickDir','out','XColor','none','YColor','none'); 
        [~,ks_stats.MeanFrTrack(this_sess)] = kstest2(remapping_raw(6).mean_rate_T1(remapping_raw(6).experiment == this_sess),remapping_raw(6).mean_rate_T2(remapping_raw(6).experiment == this_sess))
    end
 
    % peak FR
     subplot(2*num_exp,2*num_exp+2,  [num_exp+2,2*(num_exp+1)*num_exp-1]);
     [n_T1,edgesT1]= histcounts(log(remapping_raw(6).raw_peak_BAYESIAN_plfield_1),'Normalization','probability');
    n_T2= histcounts(log(remapping_raw(6).raw_peak_BAYESIAN_plfield_2),edgesT1,'Normalization','probability');
    ctrs_T1= edgesT1(1:end-1)+ 0.5*mean(diff(edgesT1));
    bar(ctrs_T1,n_T1,'FaceColor',[1 0.2 0.2],'FaceAlpha',0.5,'EdgeAlpha',0); hold on;
    bar(ctrs_T1,n_T2,'FaceColor',[0.2 0.3 1],'FaceAlpha',0.5,'EdgeAlpha',0);
    xlabel('log peak FR (Hz)'); xticks([min(edgesT1)  max(edgesT1)]);
    xticklabels(exp([min(edgesT1)  max(edgesT1)])); xtickformat('%0.1f');
    yticks([]); ylabel('probability');
    set(gca,'Box','off','TickDir','out'); 
%     legend({'Track 1','Track2'},'Box','off')
    title('Peak FR on Track')
    [~,ks_stats.PeakFrTrack_all] = kstest2(remapping_raw(6).raw_peak_BAYESIAN_plfield_1,remapping_raw(6).raw_peak_BAYESIAN_plfield_2);
    for this_sess=1:num_exp
        subplot(2*num_exp,2*num_exp+2, [2*(num_exp+1)+ (this_sess-1)*2*(num_exp+1)]);
        n_T1= histcounts(log(remapping_raw(6).raw_peak_BAYESIAN_plfield_1(remapping_raw(6).experiment == this_sess)),edgesT1,'Normalization','probability');
        n_T2= histcounts(log(remapping_raw(6).raw_peak_BAYESIAN_plfield_2(remapping_raw(6).experiment == this_sess)),edgesT1,'Normalization','probability');
        bar(ctrs_T1,n_T1,'FaceColor',[1 0.2 0.2],'FaceAlpha',0.5,'EdgeAlpha',0); hold on;
        bar(ctrs_T1,n_T2,'FaceColor',[0.2 0.3 1],'FaceAlpha',0.5,'EdgeAlpha',0);
        set(gca,'Box','off','TickDir','out','XColor','none','YColor','none'); 
       [~,ks_stats.PeakFrTrack(this_sess)] = kstest2(remapping_raw(6).raw_peak_BAYESIAN_plfield_1(remapping_raw(6).experiment == this_sess),remapping_raw(6).raw_peak_BAYESIAN_plfield_2(remapping_raw(6).experiment == this_sess))
    end

%    replay FR
    subplot(2*num_exp,2*num_exp+2,  [(num_exp)*(2*num_exp+2)+1,2*num_exp*(2*num_exp+2)-(num_exp+2)]);
    [n_T1,edgesT1]= histcounts(log(remapping_raw(6).track1_median_replay_rate),'Normalization','probability');
    n_T2= histcounts(log(remapping_raw(6).track2_median_replay_rate),edgesT1,'Normalization','probability');
    ctrs_T1= edgesT1(1:end-1)+ 0.5*mean(diff(edgesT1));
    bar(ctrs_T1,n_T1,'FaceColor',[1 0.2 0.2],'FaceAlpha',0.5,'EdgeAlpha',0); hold on;
    bar(ctrs_T1,n_T2,'FaceColor',[0.2 0.3 1],'FaceAlpha',0.5,'EdgeAlpha',0);
    xlabel('log median replay rate (Hz)'); xticks([min(edgesT1)  max(edgesT1)]);
    xticklabels(exp([min(edgesT1)  max(edgesT1)])); xtickformat('%0.1f');
    yticks([]); ylabel('probability');
    set(gca,'Box','off','TickDir','out'); 
    title('Median Replay Rate')
    [~,ks_stats.MedianFrReplay_all] = kstest2(remapping_raw(6).track1_median_replay_rate,remapping_raw(6).track2_median_replay_rate);
    for this_sess=1:num_exp
%         subplot(2*num_exp,2*num_exp+2, [5*2*(num_exp+1)+1 + num_exp+ (this_sess-1)*(2*num_exp+2)]);
        subplot(2*num_exp,2*num_exp+2, [num_exp*2*(num_exp+1)+1 + num_exp+ (this_sess-1)*(2*num_exp+2)]);        
        n_T1= histcounts(log(remapping_raw(6).track1_median_replay_rate(remapping_raw(6).experiment == this_sess)),edgesT1,'Normalization','probability');
        n_T2= histcounts(log(remapping_raw(6).track2_median_replay_rate(remapping_raw(6).experiment == this_sess)),edgesT1,'Normalization','probability');
        bar(ctrs_T1,n_T1,'FaceColor',[1 0.2 0.2],'FaceAlpha',0.5,'EdgeAlpha',0); hold on;
        bar(ctrs_T1,n_T2,'FaceColor',[0.2 0.3 1],'FaceAlpha',0.5,'EdgeAlpha',0);
        set(gca,'Box','off','TickDir','out','XColor','none','YColor','none'); 
        [~,ks_stats.MedianFrReplay(this_sess)] = kstest2(remapping_raw(6).track1_median_replay_rate(remapping_raw(6).experiment == this_sess),remapping_raw(6).track2_median_replay_rate(remapping_raw(6).experiment == this_sess))
    end
    
%    PROPORTION OF EVENTS WHERE EACH CELL IS ACTIVE
    subplot(2*num_exp,2*num_exp+2,[num_exp*(2*num_exp+2)+(num_exp+2), 2*num_exp*(2*num_exp+2)-1])
    [n_T1,edgesT1]= histcounts((remapping_raw(6).proportion_events_cell_active_track1),'Normalization','probability');
    n_T2= histcounts((remapping_raw(6).proportion_events_cell_active_track2),edgesT1,'Normalization','probability');
    ctrs_T1= edgesT1(1:end-1)+ 0.5*mean(diff(edgesT1));
    bar(ctrs_T1,n_T1,'FaceColor',[1 0.2 0.2],'FaceAlpha',0.5,'EdgeAlpha',0); hold on;
    bar(ctrs_T1,n_T2,'FaceColor',[0.2 0.3 1],'FaceAlpha',0.5,'EdgeAlpha',0);
    xlabel('proportion replay events cells are active in (Hz)'); xticks([min(edgesT1)  max(edgesT1)]);
    xticklabels(exp([min(edgesT1)  max(edgesT1)])); xtickformat('%0.1f');
    yticks([]); ylabel('probability');
    set(gca,'Box','off','TickDir','out'); 
    title('Replay Involvment')
    [~,ks_stats.ProportionEventsCellActive_all] = kstest2(remapping_raw(6).proportion_events_cell_active_track1,remapping_raw(6).proportion_events_cell_active_track2);
    for this_sess=1:num_exp
        subplot(2*num_exp,2*num_exp+2, [num_exp*2*(num_exp+1) + 2*num_exp+2+ (this_sess-1)*(2*num_exp+2)]);
        n_T1= histcounts((remapping_raw(6).proportion_events_cell_active_track1(remapping_raw(6).experiment == this_sess)),edgesT1,'Normalization','probability');
        n_T2= histcounts((remapping_raw(6).proportion_events_cell_active_track2(remapping_raw(6).experiment == this_sess)),edgesT1,'Normalization','probability');
        bar(ctrs_T1,n_T1,'FaceColor',[1 0.2 0.2],'FaceAlpha',0.5,'EdgeAlpha',0); hold on;
        bar(ctrs_T1,n_T2,'FaceColor',[0.2 0.3 1],'FaceAlpha',0.5,'EdgeAlpha',0);
        set(gca,'Box','off','TickDir','out','XColor','none','YColor','none'); 
        [~,ks_stats.ProportionEventsCellActive(this_sess)] = kstest2(remapping_raw(6).proportion_events_cell_active_track1(remapping_raw(6).experiment == this_sess),remapping_raw(6).proportion_events_cell_active_track2(remapping_raw(6).experiment == this_sess))
    end
%     
%     % REPLAY AVERAGE PEAK INSTANTEOUS FR
%     subplot(2*num_exp,2*num_exp+2,[num_exp*(2*num_exp+2)+(num_exp+2), 2*num_exp*(2*num_exp+2)-1])
%     [n_T1,edgesT1]= histcounts((remapping_raw(6).track1_mean_replay_inst_FR_nonZero),'Normalization','probability');
%     n_T2= histcounts((remapping_raw(6).track2_mean_replay_inst_FR_nonZero),edgesT1,'Normalization','probability');
%     ctrs_T1= edgesT1(1:end-1)+ 0.5*mean(diff(edgesT1));
%     bar(ctrs_T1,n_T1,'FaceColor',[1 0.2 0.2],'FaceAlpha',0.5,'EdgeAlpha',0); hold on;
%     bar(ctrs_T1,n_T2,'FaceColor',[0.2 0.3 1],'FaceAlpha',0.5,'EdgeAlpha',0);
%     xlabel('proportion replay events cells are active in (Hz)'); xticks([min(edgesT1)  max(edgesT1)]);
%     xticklabels(exp([min(edgesT1)  max(edgesT1)])); xtickformat('%0.1f');
%     yticks([]); ylabel('probability');
%     set(gca,'Box','off','TickDir','out'); 
%     title('Replay Peak Inst FR')
%     [~,ks_stats.ReplayPeakInstFR_all] = kstest2(remapping_raw(6).track1_mean_replay_inst_FR_nonZero,remapping_raw(6).track2_mean_replay_inst_FR_nonZero);
%     for this_sess=1:num_exp
%         subplot(2*num_exp,2*num_exp+2, [num_exp*2*(num_exp+1) + 2*num_exp+2+ (this_sess-1)*(2*num_exp+2)]);
%         n_T1= histcounts((remapping_raw(6).track1_mean_replay_inst_FR_nonZero(remapping_raw(6).experiment == this_sess)),edgesT1,'Normalization','probability');
%         n_T2= histcounts((remapping_raw(6).track2_mean_replay_inst_FR_nonZero(remapping_raw(6).experiment == this_sess)),edgesT1,'Normalization','probability');
%         bar(ctrs_T1,n_T1,'FaceColor',[1 0.2 0.2],'FaceAlpha',0.5,'EdgeAlpha',0); hold on;
%         bar(ctrs_T1,n_T2,'FaceColor',[0.2 0.3 1],'FaceAlpha',0.5,'EdgeAlpha',0);
%         set(gca,'Box','off','TickDir','out','XColor','none','YColor','none'); 
%         [~,ks_stats.ReplayPeakInstFR(this_sess)] = kstest2(remapping_raw(6).track1_mean_replay_inst_FR_nonZero(remapping_raw(6).experiment == this_sess),remapping_raw(6).track2_mean_replay_inst_FR_nonZero(remapping_raw(6).experiment == this_sess))
%     end
%     
%     % AUC PLACE FIELDS TACK
%       subplot(2*num_exp,2*num_exp+2,  [(num_exp)*(2*num_exp+2)+1,2*num_exp*(2*num_exp+2)-(num_exp+2)]);
%     [n_T1,edgesT1]= histcounts((remapping_raw(6).plfield_AUC_1),'Normalization','probability');
%     n_T2= histcounts((remapping_raw(6).plfield_AUC_2),edgesT1,'Normalization','probability');
%     ctrs_T1= edgesT1(1:end-1)+ 0.5*mean(diff(edgesT1));
%     bar(ctrs_T1,n_T1,'FaceColor',[1 0.2 0.2],'FaceAlpha',0.5,'EdgeAlpha',0); hold on;
%     bar(ctrs_T1,n_T2,'FaceColor',[0.2 0.3 1],'FaceAlpha',0.5,'EdgeAlpha',0);
%     xlabel('Place Field AUC'); xticks([min(edgesT1)  max(edgesT1)]);
%     xticklabels(exp([min(edgesT1)  max(edgesT1)])); xtickformat('%0.1f');
%     yticks([]); ylabel('probability');
%     set(gca,'Box','off','TickDir','out'); 
%     title('AUC')
%     [~,ks_stats.plfield_AUC_all] = kstest2(remapping_raw(6).plfield_AUC_1,remapping_raw(6).plfield_AUC_2);
%     for this_sess=1:num_exp
% %         subplot(2*num_exp,2*num_exp+2,[num_exp*(num_exp+2+ (this_sess-1)*(2*num_exp+2)]);
%         subplot(2*num_exp,2*num_exp+2, [num_exp*2*(num_exp+1)+1 + num_exp+ (this_sess-1)*(2*num_exp+2)]);
%         n_T1= histcounts((remapping_raw(6).plfield_AUC_1(remapping_raw(6).experiment == this_sess)),edgesT1,'Normalization','probability');
%         n_T2= histcounts((remapping_raw(6).plfield_AUC_2(remapping_raw(6).experiment == this_sess)),edgesT1,'Normalization','probability');
%         bar(ctrs_T1,n_T1,'FaceColor',[1 0.2 0.2],'FaceAlpha',0.5,'EdgeAlpha',0); hold on;
%         bar(ctrs_T1,n_T2,'FaceColor',[0.2 0.3 1],'FaceAlpha',0.5,'EdgeAlpha',0);
%         set(gca,'Box','off','TickDir','out','XColor','none','YColor','none'); 
%         [~,ks_stats.plfield_AUC(this_sess)] = kstest2(remapping_raw(6).plfield_AUC_1(remapping_raw(6).experiment == this_sess),remapping_raw(6).plfield_AUC_2(remapping_raw(6).experiment == this_sess))
%     end
%     
    %% number of replay events
    PRE_num_replay_T1= cellfun(@(x) size(x,2),remapping_raw(5).replay_spikes1);
    PRE_num_replay_T2= cellfun(@(x) size(x,2),remapping_raw(5).replay_spikes2);
%     PRE_num_replay_T2= PRE_num_replay_T2./PRE_num_replay_T1;
    POST_num_replay_T1= cellfun(@(x) size(x,2),remapping_raw(6).replay_spikes1);
    POST_num_replay_T2= cellfun(@(x) size(x,2),remapping_raw(6).replay_spikes2);
%     POST_num_replay_T2= POST_num_replay_T2./POST_num_replay_T1;
    RUN_num_replay_T1= cellfun(@(x) size(x,2),remapping_raw(7).replay_spikes1);
    RUN_num_replay_T2= cellfun(@(x) size(x,2),remapping_raw(7).replay_spikes2);
    
    % wilcoxon (non-normal dist)
    WCX_replay_PRE = signrank(PRE_num_replay_T1,PRE_num_replay_T2)
    WCX_replay_RUN = signrank(RUN_num_replay_T1,RUN_num_replay_T2)
    WCX_replay_POST = signrank(POST_num_replay_T1,POST_num_replay_T2)

%     rat_colours= {'r.','b.','g.','c.','m.','y.','k.'};  
%     rat_colours= {'r','b','g','c','m'};

    cmap= colormap('gray');
    cidx=linspace(1,length(cmap),length(PRE_num_replay_T2));
    rat_colours= cmap(round(cidx),:);
    
    figure('Color','w');
    ax(1)= subplot(1,3,1);
    plot([ones(size(PRE_num_replay_T1,2),1) 2*ones(size(PRE_num_replay_T2,2),1)]',[PRE_num_replay_T1' PRE_num_replay_T2']','LineWidth',1.5);
    colororder(rat_colours);
    xlim([0 3]); xticks([1 2]); xticklabels({'Track1','Track2'});
    ylabel('number of replay events'); yticks([min([PRE_num_replay_T1 PRE_num_replay_T2]) max([PRE_num_replay_T1 PRE_num_replay_T2])]);
    set(gca,'Box','off','TickDir','out'); 
    title('PRE')
    
    ax(2)=subplot(1,3,2);
    plot([ones(size(RUN_num_replay_T1,2),1) 2*ones(size(RUN_num_replay_T2,2),1)]',[RUN_num_replay_T1' RUN_num_replay_T2']','LineWidth',1.5);
    xlim([0 3]); xticks([1 2]); xticklabels({'Track1','Track2'});
    ylabel('number of replay events'); yticks([min([RUN_num_replay_T1 RUN_num_replay_T2]) max([RUN_num_replay_T1 RUN_num_replay_T2])]);
    set(gca,'Box','off','TickDir','out'); 
    title('RUN');
    colororder(rat_colours)
    
    ax(3)=subplot(1,3,3);
    plot([ones(size(POST_num_replay_T1,2),1) 2*ones(size(POST_num_replay_T2,2),1)]',[POST_num_replay_T1' POST_num_replay_T2']','LineWidth',1.5);
    xlim([0 3]); xticks([1 2]); xticklabels({'Track1','Track2'});
    ylabel('number of replay events'); yticks([min([POST_num_replay_T1 POST_num_replay_T2]) max([POST_num_replay_T1 POST_num_replay_T2])]);
    set(gca,'Box','off','TickDir','out'); 
    title('POST')
%     legend({'rat1','rat2','rat3','rat4','rat5'},'Box','off')
    linkaxes([ax(:)],'xy');
    colororder(rat_colours)
    
    
    % get replay rate
    load('folders_to_process_remapping.mat');
    master_folder= pwd;
    for this_folder= 1:size(folders,1)
        cd([master_folder '\' folders{this_folder,1}]);
        load('time_range.mat');
        time_PRE(this_folder)= diff(time_range.pre);
        time_RUN(this_folder,1)= diff(time_range.track(1).behaviour);
        time_RUN(this_folder,2)= diff(time_range.track(2).behaviour);
        time_POST(this_folder)= diff(time_range.post);
    end
    cd(master_folder);
    
    PRE_rate_replay_T1= PRE_num_replay_T1./time_PRE;
    PRE_rate_replay_T2= PRE_num_replay_T2./time_PRE;
    RUN_rate_replay_T1= RUN_num_replay_T1./time_RUN(:,1)';
    RUN_rate_replay_T2= RUN_num_replay_T2./time_RUN(:,2)';
    POST_rate_replay_T1= POST_num_replay_T1./time_POST;
    POST_rate_replay_T2= POST_num_replay_T2./time_POST;
    
    WCX_replay_rate_PRE = signrank(PRE_rate_replay_T1,PRE_rate_replay_T2)
    WCX_replay_rate_RUN = signrank(RUN_rate_replay_T1,RUN_rate_replay_T2)
    WCX_replay_rate_POST = signrank(POST_rate_replay_T1,POST_rate_replay_T2)
    
    figure('Color','w');
    ax(1)= subplot(1,3,1);
    plot([ones(size(PRE_rate_replay_T1,2),1) 2*ones(size(PRE_rate_replay_T2,2),1)]',[PRE_rate_replay_T1' PRE_rate_replay_T2']','LineWidth',1.5);
    xlim([0 3]); xticks([1 2]); xticklabels({'Track1','Track2'});
    ylabel('replay events rate (Hz)'); yticks([min([PRE_rate_replay_T1 PRE_rate_replay_T2]) max([PRE_rate_replay_T2 PRE_rate_replay_T2])]);
    set(gca,'Box','off','TickDir','out'); 
    title('PRE')
    
    ax(2)=subplot(1,3,2);
    plot([ones(size(RUN_rate_replay_T1,2),1) 2*ones(size(RUN_rate_replay_T2,2),1)]',[RUN_rate_replay_T1' RUN_rate_replay_T2']','LineWidth',1.5);
    xlim([0 3]); xticks([1 2]); xticklabels({'Track1','Track2'});
    ylabel('replay events rate (Hz)'); yticks([min([RUN_rate_replay_T1 RUN_rate_replay_T2]) max([RUN_rate_replay_T1 RUN_rate_replay_T2])]);
    set(gca,'Box','off','TickDir','out'); 
    title('RUN')
    
    ax(3)=subplot(1,3,3);
    plot([ones(size(POST_rate_replay_T1,2),1) 2*ones(size(POST_rate_replay_T2,2),1)]',[POST_rate_replay_T1' POST_rate_replay_T2']','LineWidth',1.5);
    xlim([0 3]); xticks([1 2]); xticklabels({'Track1','Track2'});
    ylabel('replay events rate (Hz)'); yticks([min([POST_rate_replay_T1 POST_rate_replay_T2]) max([POST_rate_replay_T1 POST_rate_replay_T2])]);
    set(gca,'Box','off','TickDir','out'); 
    title('POST')
    legend({'rat1','rat2','rat3','rat4','rat5'},'Box','off')
    linkaxes([ax(:)],'xy');
    
end