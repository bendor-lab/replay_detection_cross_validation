function plot_correlation_sleep_evolution(folders, method, epoch, option)

switch option
    case 'time'
        bin_edges= 60*[[0:15:60]' [30:15:90]'];
    case 'events'
        bin_edges= [[1:70:(300-70)]' [70:70:300]'];
end

folders= folders(:,1);
master_folder= 'X:\BendorLab\Drobo\Neural and Behavioural Data\Rate remapping\Data';
shuffle_folder= 'X:\BendorLab\Drobo\Neural and Behavioural Data\Rate remapping\Data\CONTROLS\Sleep Evolution' ;
temp_folders = arrayfun(@(x) [shuffle_folder '\TEMP_' folders{x}],1:length(folders),'UniformOutput',0);

for this_time_win=1:size(bin_edges,1)
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
        end
        
        load(['sorted_replay_' method '.mat']);
      switch option
                case 'time'
                % now shorten sort_replay_events to time window of interest in sleepPOST
                for this_track=1:length(sorted_replay)
                    % find ones outside window
                    evt_times_to_remove= sorted_replay(this_track).cumulative_event_time.(epoch)< bin_edges(this_time_win,1) |...
                                        sorted_replay(this_track).cumulative_event_time.(epoch)>= bin_edges(this_time_win,2);

                    % remove from struct
                    sorted_replay(this_track).cumulative_event_time.(epoch)(evt_times_to_remove)=[];
                    sorted_replay(this_track).event_time.(epoch)(evt_times_to_remove)=[];
                    sorted_replay(this_track).index.(epoch)(evt_times_to_remove)=[];
                    sorted_replay(this_track).ref_index.(epoch)(evt_times_to_remove)=[];
                end
                
          case 'events'
              % find 20 events T1+T2
              [~,all_events_idx]= sort([sorted_replay(1).event_time.(epoch) sorted_replay(2).event_time.(epoch)]);
              track_idx= [ones(1,size(sorted_replay(1).event_time.(epoch),2)) 2*ones(1,size(sorted_replay(1).event_time.(epoch),2))];
              if bin_edges(this_time_win,2) >= length(all_events_idx) && bin_edges(this_time_win,1) >= length(all_events_idx) 
                  idx_to_keep=[];
              elseif bin_edges(this_time_win,2) > length(all_events_idx) && bin_edges(this_time_win,1) < length(all_events_idx) 
                  idx_to_keep= all_events_idx(bin_edges(this_time_win,1):end);
              else
                idx_to_keep= all_events_idx(bin_edges(this_time_win,1):bin_edges(this_time_win,2));
              end
              evt_times_to_keep_T1= find(ismember([1:length(track_idx)],idx_to_keep) & track_idx==1);
              evt_times_to_remove_T2= find(ismember([1:length(track_idx)],idx_to_keep) & track_idx==2) - find(track_idx==1,1,'last');
              %T1
              sorted_replay(1).cumulative_event_time.(epoch)= sorted_replay(1).cumulative_event_time.(epoch)(evt_times_to_keep_T1);
              sorted_replay(1).event_time.(epoch)= sorted_replay(1).event_time.(epoch)(evt_times_to_keep_T1);
              sorted_replay(1).index.(epoch)= sorted_replay(1).index.(epoch)(evt_times_to_keep_T1);
              sorted_replay(1).ref_index.(epoch)= sorted_replay(1).ref_index.(epoch)(evt_times_to_keep_T1);
              %T2
              sorted_replay(2).cumulative_event_time.(epoch)= sorted_replay(2).cumulative_event_time.(epoch)(evt_times_to_remove_T2);
              sorted_replay(2).event_time.(epoch)= sorted_replay(2).event_time.(epoch)(evt_times_to_remove_T2);
              sorted_replay(2).index.(epoch)= sorted_replay(2).index.(epoch)(evt_times_to_remove_T2);
              sorted_replay(2).ref_index.(epoch)= sorted_replay(2).ref_index.(epoch)(evt_times_to_remove_T2);
      end
        if ~isempty(sorted_replay(1).cumulative_event_time.(epoch)) || ~isempty(sorted_replay(2).cumulative_event_time.(epoch))
            num_sess_this_win(this_folder,this_time_win)= 1;
        else
            num_sess_this_win(this_folder,this_time_win)= 0;
        end
        save(['sorted_replay_' method '.mat'],'sorted_replay');
        
        
    end
    cd ..
    [remapping, remapping_raw] = rate_remapping_TRACK_PAIRS(temp_folders,'wcorr',0);
    cd(master_folder)
    [pval(this_time_win),Fstat(this_time_win), s_intcp(this_time_win,:)] = plot_rate_remapping_NEW('use_mat',remapping,'x_var',{'place_field_diff'},'y_var',....
    {'mean_max_FR_replay_diff'},'epochs',{epoch},'subset','stable cells laps');
    f=gcf;
    ax=f.Children; hs = findobj(ax,'Type','Scatter');
    num_data_points(this_time_win)= numel(hs.XData);
    close gcf;
    epoch_idx= strcmp([remapping.epoch],epoch);
    num_events_T1_this_win(this_time_win)= sum(arrayfun(@(x) size(remapping_raw(epoch_idx).replay_spikes1{x},2),1:length(remapping_raw(epoch_idx).replay_spikes1)));
    num_events_T2_this_win(this_time_win)= sum(arrayfun(@(x) size(remapping_raw(epoch_idx).replay_spikes2{x},2),1:length(remapping_raw(epoch_idx).replay_spikes2)));

end
cd(shuffle_folder)
sleep_evolution_corr.num_sess_this_win= num_sess_this_win;
sleep_evolution_corr.num_data_points= num_data_points;
% sleep_evolution_corr.num_events_T1_this_win= num_events_T1_this_win;
% sleep_evolution_corr.num_events_T2_this_win= num_events_T2_this_win;
sleep_evolution_corr.pval= pval;
sleep_evolution_corr.Fstat= Fstat;
sleep_evolution_corr.slope= s_intcp(:,1);
sleep_evolution_corr.intercept= s_intcp(:,2);
sleep_evolution_corr.bin_edges= bin_edges;
sleep_evolution_corr.method= option;
save(['sleep_evolution_' epoch '_corr.mat'],'sleep_evolution_corr');

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
[~,~,bin_id]= histcounts(sleep_evolution_corr.pval,new_steps);

sleep_time_edges= sleep_evolution_corr.bin_edges;
sleep_ctrs= sleep_time_edges(1:end,1)+0.5*mean(diff(sleep_time_edges(:,2)));

f3= figure('Color','w','Name',epoch);
p1= plot(sleep_ctrs,sleep_evolution_corr.Fstat,'k','LineWidth',1.5);
hold on;
scatter(sleep_ctrs,sleep_evolution_corr.Fstat,50,colmap(bin_id,:),'filled');
colormap(colmap);
%ylim([-30 40]); yticks([0:10:40])
ylim([0 60]); yticks([0:10:60])
ylabel('F-Statistic')
ytickformat('%i')
yyaxis right
% tot_events= (sleep_evolution_corr.num_events_T1_this_win + sleep_evolution_corr.num_events_T2_this_win);
% p2 = plot(sleep_ctrs,100*sleep_evolution_corr.num_events_T1_this_win./tot_events,'Color', [0.7 0.7 0.7],'LineWidth',1.5);
% sum_events= tot_events./sum(tot_events);
% p3 = plot(sleep_ctrs,100*sum_events,'Color', [0.7 0.7 0.7],'LineWidth',2,'LineStyle','-');
p3 = plot(sleep_ctrs,sleep_evolution_corr.num_data_points,'Color', [0.7 0.7 0.7],'LineWidth',2,'LineStyle','-');
p2= plot(sleep_ctrs,sum(sleep_evolution_corr.num_sess_this_win),'Color', [0.7 0.7 0.7],'LineWidth',2,'LineStyle','-');
 set(gca,{'ycolor'},{[0.4 0.4 0.4]})
ylabel('Number of Cells','Rotation',270);  %ylim([0 15]); yticks([0 5]); yticklabels({0,'max'})
% yticks(round(100*[min(sum_events) max(sum_events)]) ) %...
%     min(sleep_evolution_corr.num_events_T1_this_win./tot_events) max(sleep_evolution_corr.num_events_T1_this_win./tot_events)])); 
% ylim([0 30]);
xticks(round(linspace(sleep_ctrs(1,1),sleep_ctrs(end,1),length(sleep_ctrs))));
% xticklabels(xticks./60);
% xlim([sleep_time_edges(1,1) sleep_ctrs(end)+200])
xlim([sleep_time_edges(1,1) sleep_ctrs(end)+20]);
% xlabel('sleep time (min)'); 
xlabel('Event Block (each session)');
% legend([p1,p3],{'FStat','% total events'},'Box','off','Location','northeast');
run_format_settings(f3) 
% legend([p1,p2,p3],{'FStat','number of events T1/(T1+T2)','% total events'},'Box','off','Location','northeast');
legend([p1,p2,p3],{'FStat','Number of Sessions','Number of Cells'},'Box','off','Location','northeast');
h = colorbar;
colormap(colmap)
run_format_settings(f3,'match_ax')
h.Ticks= [0 find(new_steps==0.05)/size(colmap,1) 1]; 
h.TickLabels = [0 0.05 1];
ylabel(h,'pval','Rotation',270);
%
    
% 
% yyaxis right
% plot(sleep_ctrs,sleep_evolution_corr.pval,'r-','LineWidth',1.5,'MarkerSize',20);
% plot(sleep_ctrs,sleep_evolution_corr.slope.*sum(gausswin(0.25/(1/1000))),'Color',[0.5 0.5 0.5],'LineWidth',1.5,'MarkerSize',20,'Marker','.');
% ylim([-0.01  0.12]); 
% legend({'FStat','p value','slope'},'Box','off','Location','southwest')
% % set(gca,'YColor','r')
% ylabel('pvalue, slope')


end