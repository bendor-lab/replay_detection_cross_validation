function plot_correlation_RUN_evolution(folders, method, option)
%methods: wcorr or spearman

switch option
    case 'time'
        bin_edges= 60*[[0:1:24-1]' [1:1:24]'];
        num_bins = size(bin_edges,1);
    case {'laps','cumulative_laps'}
        num_laps = 20;
        jump_laps = 1;
        num_bins = floor(20/jump_laps);
    case 'events'
        bin_edges= [[1:10:91]' [10:10:100]'];
        num_bins = size(bin_edges,1);
end

if isempty(folders)
    load('X:\BendorLab\Drobo\Neural and Behavioural Data\Rate remapping\Data\folders_to_process_remapping.mat')
    folders= folders(:,1);
else
    folders= folders(:,1);
end
master_folder= 'X:\BendorLab\Drobo\Neural and Behavioural Data\Rate remapping\Data';
shuffle_folder= 'X:\BendorLab\Drobo\Neural and Behavioural Data\Rate remapping\Data\CONTROLS\RUN Evolution' ;
temp_folders = arrayfun(@(x) [shuffle_folder '\TEMP_' folders{x}],1:length(folders),'UniformOutput',0);

for this_time_win = 1 : num_bins
    for this_folder=1:length(folders)
        if exist(temp_folders{this_folder})~=7
            mkdir(temp_folders{this_folder})
        end
        cd(temp_folders{this_folder});
        copyfile(['..\..\..\' folders{this_folder} '\sorted_replay_' method '.mat']);
        if ~exist('significant_replay_events_wcorr.mat') % no need to copy everything again
            if exist(['..\..\..\' folders{this_folder} '\significant_replay_events_wcorr.mat'])==2
                copyfile(['..\..\..\' folders{this_folder} '\significant_replay_events_wcorr.mat']);
            elseif exist(['..\..\..\' folders{this_folder} '\significant_replay_events_wcorr_individual_exposures.mat'])==2
                 copyfile(['..\..\..\' folders{this_folder} '\significant_replay_events_wcorr_individual_exposures.mat']);
            end
            copyfile(['..\..\..\' folders{this_folder} '\extracted_place_fields_BAYESIAN.mat']);
            copyfile(['..\..\..\' folders{this_folder} '\lap_times.mat']);
        end
        
        load(['sorted_replay_' method '.mat']);
        switch option
            case 'laps'
                load('lap_times.mat')
                for this_track=1:length(sorted_replay)
                    bin_edges{this_track} = [lap_times(this_track).start(1:jump_laps:num_laps-(jump_laps-1))' lap_times(this_track).end(jump_laps:jump_laps:num_laps)'];
                    % find ones outside window
                    evt_times_to_remove= sorted_replay(this_track).event_time.track(this_track).behaviour < bin_edges{this_track}(this_time_win,1) |...
                        sorted_replay(this_track).event_time.track(this_track).behaviour>= bin_edges{this_track}(this_time_win,2);
                    
                    % remove from struct
                    sorted_replay(this_track).cumulative_event_time.track(this_track).behaviour(evt_times_to_remove)=[];
                    sorted_replay(this_track).event_time.track(this_track).behaviour(evt_times_to_remove)=[];
                    sorted_replay(this_track).index.track(this_track).behaviour(evt_times_to_remove)=[];
                    sorted_replay(this_track).ref_index.track(this_track).behaviour(evt_times_to_remove)=[];
                end
                
            case 'cumulative_laps'
                load('lap_times.mat')
                for this_track=1:length(sorted_replay)
                    bin_edges{this_track} = [ones(num_laps,1)*lap_times(this_track).start(1) lap_times(this_track).end(jump_laps:jump_laps:num_laps)'];
                    % find ones outside window
                    evt_times_to_remove= sorted_replay(this_track).event_time.track(this_track).behaviour < bin_edges{this_track}(this_time_win,1) |...
                        sorted_replay(this_track).event_time.track(this_track).behaviour>= bin_edges{this_track}(this_time_win,2);
                    
                    % remove from struct
                    sorted_replay(this_track).cumulative_event_time.track(this_track).behaviour(evt_times_to_remove)=[];
                    sorted_replay(this_track).event_time.track(this_track).behaviour(evt_times_to_remove)=[];
                    sorted_replay(this_track).index.track(this_track).behaviour(evt_times_to_remove)=[];
                    sorted_replay(this_track).ref_index.track(this_track).behaviour(evt_times_to_remove)=[];
                end
                
            case 'time'
                load('lap_times.mat')
                cumulative_time{1} =compute_cumulative_time([lap_times(1).start; lap_times(1).end]);
                cumulative_time{2} =compute_cumulative_time([lap_times(2).start; lap_times(2).end]);
                % now shorten sort_replay_events to time window of interest in sleepPOST
                for this_track=1:2%length(sorted_replay)
                    % find ones outside window
                    evt_times_to_remove= sorted_replay(this_track).cumulative_event_time.track(this_track).behaviour < bin_edges(this_time_win,1) |...
                        sorted_replay(this_track).cumulative_event_time.track(this_track).behaviour>= bin_edges(this_time_win,2);
                    % Find lap IDs during this time window
                    [~,closest_lap_end] = min(abs(cumulative_time{1,this_track}(2,:) - bin_edges(this_time_win,2)));
                    if this_time_win == 1
                        if this_track == 1
                            T1_timeWin_LapID(this_folder,this_time_win) = closest_lap_end;
                            T1_timeWin_num_laps(this_folder,this_time_win) = closest_lap_end;
                        else
                            T2_timeWin_LapID(this_folder,this_time_win) = closest_lap_end;
                            T2_timeWin_num_laps(this_folder,this_time_win) = closest_lap_end;
                        end
                    else
                        if this_track == 1
                            T1_timeWin_LapID(this_folder,this_time_win) = closest_lap_end;
                            T1_timeWin_num_laps(this_folder,this_time_win) = closest_lap_end - T1_timeWin_num_laps(this_folder,this_time_win-1);
                        else
                            T2_timeWin_LapID(this_folder,this_time_win) = closest_lap_end;
                            T2_timeWin_num_laps(this_folder,this_time_win) = closest_lap_end - T2_timeWin_num_laps(this_folder,this_time_win-1);
                        end
                    end
                    % remove from struct
                    sorted_replay(this_track).cumulative_event_time.track(this_track).behaviour(evt_times_to_remove)=[];
                    sorted_replay(this_track).event_time.track(this_track).behaviour(evt_times_to_remove)=[];
                    sorted_replay(this_track).index.track(this_track).behaviour(evt_times_to_remove)=[];
                    sorted_replay(this_track).ref_index.track(this_track).behaviour(evt_times_to_remove)=[];
                end
                
            case 'events'
                load('lap_times.mat')
                cumulative_time{1} =compute_cumulative_time([lap_times(1).start; lap_times(1).end]);
                cumulative_time{2} =compute_cumulative_time([lap_times(2).start; lap_times(2).end]);
                % find 20 events for each track separate
                for t = 1 : 2%length(sorted_replay)
                    local_events_idx = 1: length(sorted_replay(t).event_time.track(t).behaviour);
                    %track_idx= [ones(1,size(sorted_replay(1).event_time.(epoch),2)) 2*ones(1,size(sorted_replay(1).event_time.(epoch),2))];
                    if bin_edges(this_time_win,2) >= length(local_events_idx) && bin_edges(this_time_win,1) >= length(local_events_idx)
                        idx_to_keep=[];
                    elseif bin_edges(this_time_win,2) > length(local_events_idx) && bin_edges(this_time_win,1) < length(local_events_idx)
                        idx_to_keep= local_events_idx(bin_edges(this_time_win,1):end);
                    else
                        idx_to_keep= local_events_idx(bin_edges(this_time_win,1):bin_edges(this_time_win,2));
                    end
                    sorted_replay(t).cumulative_event_time.track(t).behaviour = sorted_replay(t).cumulative_event_time.track(t).behaviour(idx_to_keep);
                    sorted_replay(t).event_time.track(t).behaviour = sorted_replay(t).event_time.track(t).behaviour(idx_to_keep);
                    sorted_replay(t).index.track(t).behaviour = sorted_replay(t).index.track(t).behaviour(idx_to_keep);
                    sorted_replay(t).ref_index.track(t).behaviour = sorted_replay(t).ref_index.track(t).behaviour(idx_to_keep);
                    % Find lap IDs during this time window
                    if ~isempty(sorted_replay(t).cumulative_event_time.track(t).behaviour)
                        [~,closest_lap_end] = min(abs(cumulative_time{1,t}(2,:) - sorted_replay(t).cumulative_event_time.track(t).behaviour(end)));
                        if this_time_win == 1
                            if t == 1
                                T1_timeWin_LapID(this_folder,this_time_win) = closest_lap_end;
                                T1_timeWin_num_laps(this_folder,this_time_win) = closest_lap_end;
                            else
                                T2_timeWin_LapID(this_folder,this_time_win) = closest_lap_end;
                                T2_timeWin_num_laps(this_folder,this_time_win) = closest_lap_end;
                            end
                        else
                            if t == 1
                                T1_timeWin_LapID(this_folder,this_time_win) = closest_lap_end;
                                T1_timeWin_num_laps(this_folder,this_time_win) = closest_lap_end - T1_timeWin_num_laps(this_folder,this_time_win-1);
                            else
                                T2_timeWin_LapID(this_folder,this_time_win) = closest_lap_end;
                                T2_timeWin_num_laps(this_folder,this_time_win) = closest_lap_end - T2_timeWin_num_laps(this_folder,this_time_win-1);
                            end
                        end
                    end
                end
      end
      for t = 1 : length(sorted_replay)
          if ~isempty(sorted_replay(t).cumulative_event_time.track(t).behaviour)
              num_sess_this_win(this_folder,this_time_win)= 1;
          else
              num_sess_this_win(this_folder,this_time_win)= 0;
          end
      end
      save(['sorted_replay_' method '.mat'],'sorted_replay');
      
      
    end
    cd ..
    [remapping, remapping_raw] = rate_remapping_TRACK_PAIRS(temp_folders,'wcorr',0);
    cd(master_folder) 
    epoch_idx= strcmp([remapping.epoch],'RUN');
    if ~isempty(remapping(epoch_idx).place_field_diff)
        [pval(this_time_win),Fstat(this_time_win), ~] = plot_rate_remapping_NEW('use_mat',remapping,'x_var',{'place_field_diff'},'y_var',....
            {'mean_max_FR_replay_diff'},'epochs',{'RUN'},'subset','stable cells laps');%s_intcp(this_time_win,:)
        num_events_T1_this_win(this_time_win)= sum(arrayfun(@(x) size(remapping_raw(epoch_idx).replay_spikes1{x},2),1:length(remapping_raw(epoch_idx).replay_spikes1)));
        num_events_T2_this_win(this_time_win)= sum(arrayfun(@(x) size(remapping_raw(epoch_idx).replay_spikes2{x},2),1:length(remapping_raw(epoch_idx).replay_spikes2)));
        f=gcf;
        ax=f.Children; hs = findobj(ax,'Type','Scatter');
        num_data_points(this_time_win)= numel(hs.XData);
        close gcf;
    else
        pval(this_time_win) = NaN;
        Fstat(this_time_win) = NaN;
        s_intcp(this_time_win,:) = NaN;
    end
    

end
cd(shuffle_folder)
run_evolution_corr.num_sess_this_win= num_sess_this_win;
run_evolution_corr.num_events_T1_this_win= num_events_T1_this_win;
run_evolution_corr.num_events_T2_this_win= num_events_T2_this_win;
run_evolution_corr.num_data_points= num_data_points;
run_evolution_corr.pval= pval;
run_evolution_corr.Fstat= Fstat;
%run_evolution_corr.slope= s_intcp(:,1);
%run_evolution_corr.intercept= s_intcp(:,2);
run_evolution_corr.bin_edges= bin_edges;
run_evolution_corr.method= option;
if exist('T1_timeWin_LapID','var')
    run_evolution_corr.T1_timeWin_LapID= T1_timeWin_LapID;
    run_evolution_corr.T1_timeWin_num_laps= T1_timeWin_num_laps;
    run_evolution_corr.T2_timeWin_LapID= T2_timeWin_LapID;
    run_evolution_corr.T2_timeWin_num_laps= T2_timeWin_num_laps;
end
save(['run_evolution_' option '_corr.mat'],'run_evolution_corr');

% Delete temporary folders
arrayfun(@(x) rmdir(temp_folders{x},'s'),1:length(temp_folders));

% make fig
% create colormap
steps = 0:0.001:0.05;
colmap = zeros(length(steps),3);
colmap(:,1) = linspace(1,0.7,length(steps)); 
colmap(:,2) = linspace(0,0,length(steps));
colmap(:,3) = linspace(0,0.7,length(steps));
added_steps= linspace(0.06,1,5);
colmap2 = zeros(length(added_steps),3);
colmap2(:,1)= linspace(0.5,0,length(added_steps));
colmap2(:,2) = linspace(0,0,length(added_steps));
colmap2(:,3) = linspace(0.5,1,length(added_steps));

colmap= ([colmap ; colmap2]);
new_steps= [steps added_steps];
[~,~,bin_id]= histcounts(run_evolution_corr.pval,new_steps);
if any(ismember(bin_id,0))
    pval_colmap = colmap(bin_id(~ismember(bin_id,0)),:);
    nan_idx = find(ismember(bin_id,0));
    if nan_idx == 1
        pval_colmap = [0.3 0.3 0.3; pval_colmap];
    elseif nan_idx == length(pval)
        pval_colmap = [pval_colmap; 0.3 0.3 0.3];
    else
        pre_nan_idcs = 1: nan_idx-1;
        post_nan_idcs =  nan_idx+1:length(pval);
        pval_colmap = [pval_colmap(pre_nan_idcs,:); 0.3 0.3 0.3;pval_colmap(post_nan_idcs,:)];
    end    
else
    pval_colmap = colmap(bin_id,:);
end
    
if strcmp(option,'laps') | strcmp(option,'cumulative_laps')
    run_time_edges = 1:jump_laps:num_laps;
    run_ctrs = run_time_edges;
else
    run_time_edges= run_evolution_corr.bin_edges;
    run_ctrs= run_time_edges(1:end,1)+0.5*mean(diff(run_time_edges(:,2)));
end

f3= figure('Color','w','Name','RUN');
p1= plot(run_ctrs,run_evolution_corr.Fstat,'k','LineWidth',1.5);
hold on;
scatter(run_ctrs,run_evolution_corr.Fstat,50,pval_colmap,'filled');
colormap(colmap);
ylabel('F Statistic')
ytickformat('%i')

yyaxis right
p2 = plot(run_ctrs,sum(run_evolution_corr.num_sess_this_win),'Color', [0.7 0.7 0.7],'LineWidth',2,'LineStyle','-');
set(gca,{'ycolor'},{[0.4 0.4 0.4]})
ylabel('number of sessions','Rotation',270); ylim([0 100]); 

tot_events= (run_evolution_corr.num_events_T1_this_win + run_evolution_corr.num_events_T2_this_win);
% p2 = plot(run_ctrs,100*run_evolution_corr.num_events_T1_this_win./tot_events,'Color', [0.7 0.7 0.7],'LineWidth',1.5);
sum_events= tot_events./sum(tot_events);
p3 = plot(run_ctrs,100*sum_events,'Color', [0.7 0.5 0.7],'LineWidth',2,'LineStyle','-');
p4 =  plot(run_ctrs,run_evolution_corr.num_data_points,'Color', [0.5 0.7 0.7],'LineWidth',2,'LineStyle','-');


if strcmp(option,'events')
    xlabel('event number');
    xticks(round(linspace(run_ctrs(1,1),run_ctrs(end,1),3)));
    xlim([run_time_edges(1,1) run_ctrs(end)+10]);
elseif strcmp(option,'time')
    xlabel('Time')
    xticks(round(linspace(run_ctrs(1,1),run_ctrs(end,1),3)));
    xlim([run_time_edges(1,1) run_ctrs(end)+10]);
else
    xlabel('Laps')
    xticks(0:5:num_laps);
    xticklabels(0:5:num_laps);
end

legend([p1,p2,p4,p3],{'FStat','number of sessions','Number of data points','Proportion of total events'},'Box','off','Location','northeast');
h = colorbar;
colormap(colmap)
run_format_settings(f3,'match_ax')
h.Ticks= [0 find(new_steps==0.05)/size(colmap,1) 1]; 
h.TickLabels = [0 0.05 1];
ylabel(h,'pval','Rotation',270);
yyaxis right
yticks([min(ylim):20: max(ylim)]); 

if  ~strcmp(option,'laps')
    figure
    boxplot([run_evolution_corr.T1_timeWin_LapID; run_evolution_corr.T2_timeWin_LapID],'PlotStyle','traditional','Color',[.3 .3 .3],'LabelOrientation','horizontal','Widths',0.5);%,'BoxStyle','filled'
    a = get(get(gca,'children'),'children');   % Get the handles of all the objects
    tt = get(a,'tag');   % List the names of all the objects
    idx = find(strcmpi(tt,'box')==1);  % Find Box objects
    boxes =a(idx);
    set(boxes,'LineWidth',1); % Set width
    box off
    hold on
    for jj = 1 : size(run_evolution_corr.T1_timeWin_LapID,2)
        plot(jj,run_evolution_corr.T1_timeWin_LapID(:,jj),'o','MarkerFaceColor','w','MarkerEdgeColor','b','LineWidth',1.5)
        plot(jj,run_evolution_corr.T2_timeWin_LapID(:,jj),'o','MarkerFaceColor','w','MarkerEdgeColor','r','LineWidth',1.5)
    end
    set(gcf,'Color','w')
    set(gca,'FontSize',20)
    xlabel('Time bin number')
    ylabel('Lap ID')
end


end


function cumulative_time=compute_cumulative_time(time)
% Takes start and end times for a specific period (e.g. PRE sleep), normalized it to the start of the sleep epoch within the sleep period and 
% then calculates cummulative sum for the start and end timestamps (stacks all epochs together)

cumulative_time=time-time(1,:); %normalize by epoch start
cumulative_time(2,:)=cumsum(cumulative_time(2,:));%cummulative sum of stop times
cumulative_time(1,2:end)=cumulative_time(2,1:(end-1));
end