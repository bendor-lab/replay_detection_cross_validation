function replay_firing_stats= get_firing_stats_replay(folders,method)
% calculates firing distributions during replay for cells included in
% rate_remapping analysis for each epoch
% if epochs are added in remapping struct, need to update this code as well
% distributions returned: time between first and last spike in event (cell-wise)
%                                       ISI if more than one spike
%                                       number of spikes
%                                       firing rate


load(['rate_remapping_analysis_TRACK_PAIRS_' method '.mat']);


tracks_to_compare= [1 2];
num_spikes_edges= [0:1:50];
num_spikes_ctrs= num_spikes_edges(1:end-1)+0.5*mean(diff(num_spikes_edges));
FR_edges= [0:1:100];
FR_ctrs= FR_edges(1:end-1)+0.5*mean(diff(FR_edges));
dur_spiking_edges= [50:5:500]*0.001;
dur_spiking_ctrs= dur_spiking_edges(1:end-1)+0.5*mean(diff(dur_spiking_edges));
ISI_edges= [0.003:0.001:0.500]; % 3ms to 500ms...
ISI_ctrs= ISI_edges(1:end-1)+0.5*mean(diff(ISI_edges));

replay_firing_stats= struct;
for i=1:size(remapping,1)
    replay_firing_stats(i).experiment= remapping(i).experiment;
    replay_firing_stats(i).folder= remapping(i).folder;
    replay_firing_stats(i).epoch= remapping(i).epoch;
    replay_firing_stats(i).ID_active_cells_during_replay= remapping(i).ID_active_cells_during_replay;
    replay_firing_stats(i).new_ID= remapping(i).new_ID;
end



master_folder= pwd;
for this_folder= 1:length(folders)
    cd([master_folder '\' folders{this_folder}]);
    
    sorted_replay = sort_replay_events([],'wcorr')
    significant_replay_events = number_of_significant_replays(0.05,3, 'wcorr', [])
    
    for epoch = 1:6  %sleep PRE &  POST, awake PRE & POST
        if epoch==1
            events1 = significant_replay_events.track(tracks_to_compare(1)); %info from events in T1
%             index1  = intersect(sorted_replay(tracks_to_compare(1)).index.sleepPRE,events1.ref_index); %indices replay events PRE for T1
            index1  = sorted_replay(tracks_to_compare(1)).index.sleepPRE; %indices replay events PRE for T1
            events2 = significant_replay_events.track(tracks_to_compare(2)); %info from events in T2
            index2  = sorted_replay(tracks_to_compare(2)).index.sleepPRE; %indices replay events PRE for T2
        elseif epoch==2
            index1  = sorted_replay(tracks_to_compare(1)).index.sleepPOST; %indices replay events POST for T1
            events1 = significant_replay_events.track(tracks_to_compare(1));
            index2  = sorted_replay(tracks_to_compare(2)).index.sleepPOST; %indices replay events POST for T2
            events2 = significant_replay_events.track(tracks_to_compare(2));
        elseif epoch==3
            index1  = sorted_replay(tracks_to_compare(1)).index.awakePRE; 
            events1 = significant_replay_events.track(tracks_to_compare(1)); 
            index2  = sorted_replay(tracks_to_compare(2)).index.awakePRE; 
            events2 = significant_replay_events.track(tracks_to_compare(2)); 
        elseif epoch==4
            index1  = sorted_replay(tracks_to_compare(1)).index.awakePOST; 
            events1 = significant_replay_events.track(tracks_to_compare(1)); 
            index2  = sorted_replay(tracks_to_compare(2)).index.awakePOST; 
            events2 = significant_replay_events.track(tracks_to_compare(2)); 
        elseif epoch==5
            index1  = sorted_replay(tracks_to_compare(1)).index.PRE; 
            events1 = significant_replay_events.track(tracks_to_compare(1)); 
            index2  = sorted_replay(tracks_to_compare(2)).index.PRE; 
            events2 = significant_replay_events.track(tracks_to_compare(2)); 
        elseif epoch==6
            index1  = sorted_replay(tracks_to_compare(1)).index.POST; 
            events1 = significant_replay_events.track(tracks_to_compare(1)); 
            index2  = sorted_replay(tracks_to_compare(2)).index.POST; 
            events2 = significant_replay_events.track(tracks_to_compare(2)); 
        end
        
        % now for each cell get distributions (only cells active during replay)
        session_indices= remapping(epoch).experiment == find(strcmp(remapping(epoch).folder,folders{this_folder}));
        cell_indices= remapping(epoch).new_ID(session_indices);
        for this_cell= 1:length(cell_indices)
            cell_idx= find(remapping(epoch).new_ID ==  cell_indices(this_cell));
            cell_id= remapping(epoch).ID_active_cells_during_replay(cell_idx);
            
            dur_spiking= []; num_spikes= []; ISI_spikes= []; FR= [];
            % for each replay event where cell is active, and each track
            for j=1:length(index1)
                % find spikes
                spikes_tmp= find(events1.spikes{index1(j)}(:,1)== cell_id);
                spike_ts_tmp= events1.spikes{index1(j)}(spikes_tmp,2);
                num_spikes=  [num_spikes; [length(spikes_tmp) 1]];
                
                FR= [FR;  [length(spikes_tmp)./events1.event_duration(index1(j))  1]];
                if length(spike_ts_tmp)>1 
                    % duration last - first spike
                    dur_spiking= [dur_spiking; [spike_ts_tmp(end) - spike_ts_tmp(1) 1]];
                    
                    % ISI
                    ISI_spikes= [ISI_spikes; [diff(spike_ts_tmp) ones(length(diff(spike_ts_tmp)),1)]]; 
                end
            end
            for j=1:length(index2)
                % find spikes
                spikes_tmp= find(events2.spikes{index2(j)}(:,1)== cell_id);
                spike_ts_tmp= events2.spikes{index2(j)}(spikes_tmp,2);
                
                num_spikes=  [num_spikes; [length(spikes_tmp) 2]];
                FR= [FR;  [length(spikes_tmp)./events2.event_duration(index2(j))  2]];
                if length(spike_ts_tmp)>1
                    % duration last - first spike
                    dur_spiking= [dur_spiking; [spike_ts_tmp(end) - spike_ts_tmp(1) 2]];
                    
                    % ISI
                    ISI_spikes= [ISI_spikes; [diff(spike_ts_tmp) 2*ones(length(diff(spike_ts_tmp)),1)]]; 
                end
            end
            % now get cell-specific distributions
            % first column is T1 events, second is T2 events, third is all
            % events
            replay_firing_stats(epoch).num_spikes{cell_idx,1}= histcounts(num_spikes(num_spikes(:,2)==1,1),num_spikes_edges,'Normalization','probability');
            replay_firing_stats(epoch).num_spikes{cell_idx,2}= histcounts(num_spikes(num_spikes(:,2)==2,1),num_spikes_edges,'Normalization','probability');
            replay_firing_stats(epoch).num_spikes{cell_idx,3}= histcounts(num_spikes(:,1),num_spikes_edges,'Normalization','probability');
            replay_firing_stats(epoch).FR{cell_idx,1}= histcounts(FR(FR(:,2)==1,1),FR_edges,'Normalization','probability');
            replay_firing_stats(epoch).FR{cell_idx,2}= histcounts(FR(FR(:,2)==2,1),FR_edges,'Normalization','probability');
            replay_firing_stats(epoch).FR{cell_idx,3}= histcounts(FR(:,1),FR_edges,'Normalization','probability');
            if ~isempty(dur_spiking)
                replay_firing_stats(epoch).max_spike_spacing{cell_idx,1}= histcounts(dur_spiking(dur_spiking(:,2)==1,1),dur_spiking_edges,'Normalization','probability');
                replay_firing_stats(epoch).max_spike_spacing{cell_idx,2}= histcounts(dur_spiking(dur_spiking(:,2)==2,1),dur_spiking_edges,'Normalization','probability');
                replay_firing_stats(epoch).max_spike_spacing{cell_idx,3}= histcounts(dur_spiking(:,1),dur_spiking_edges,'Normalization','probability');
                replay_firing_stats(epoch).ISI{cell_idx,1}= histcounts(ISI_spikes(ISI_spikes(:,2)==1,1),ISI_edges,'Normalization','probability');
                replay_firing_stats(epoch).ISI{cell_idx,2}= histcounts(ISI_spikes(ISI_spikes(:,2)==2,1),ISI_edges,'Normalization','probability');
                replay_firing_stats(epoch).ISI{cell_idx,3}= histcounts(ISI_spikes(:,1),ISI_edges,'Normalization','probability');
            else
                for jj=1:3
                    replay_firing_stats(epoch).max_spike_spacing{cell_idx,jj}= NaN(size(dur_spiking_ctrs));
                    replay_firing_stats(epoch).ISI{cell_idx,jj}= NaN(size(ISI_ctrs));
                end
            end
        end

    end

end
    for epoch = 1:6
        for jj=1:3
            % get average over cells
             replay_firing_stats(epoch).num_spikes_all_cells{jj}= nanmean(cell2mat(vertcat(replay_firing_stats(epoch).num_spikes(:,jj)))); 
             replay_firing_stats(epoch).num_spikes_ctrs= num_spikes_ctrs;
             replay_firing_stats(epoch).FR_all_cells{jj}= nanmean(cell2mat(vertcat(replay_firing_stats(epoch).FR(:,jj)))); 
             replay_firing_stats(epoch).FR_ctrs= FR_ctrs;
             replay_firing_stats(epoch).max_spike_spacing_all_cells{jj}= nanmean(cell2mat(vertcat(replay_firing_stats(epoch).max_spike_spacing(:,jj)))); 
             replay_firing_stats(epoch).max_spike_spacing_ctrs= dur_spiking_ctrs;
             replay_firing_stats(epoch).ISI_all_cells{jj}= nanmean(cell2mat(vertcat(replay_firing_stats(epoch).ISI(:,jj)))); 
             replay_firing_stats(epoch).ISI_ctrs= ISI_ctrs;
        end
    end
    cd(master_folder);
    
    save('.\Tables\replay_firing_stats.mat','replay_firing_stats');
    
    % number of consecutive bins with spikes?
    
%% figures
title_label= {'Track1 events','Track2 events','All Events'};
figure('color','w'); % POST
for jj=1:3
    subplot(1,3,jj);
    bar(replay_firing_stats(2).max_spike_spacing_ctrs,replay_firing_stats(2).max_spike_spacing_all_cells{jj});
    xlabel('time between first and last spikes, averaged all cells (s)')
    ylabel('probability')
    title(title_label{jj});
end

figure('color','w'); % POST
for jj=1:3
    subplot(1,3,jj);
    bar(replay_firing_stats(2).ISI_ctrs,replay_firing_stats(2).ISI_all_cells{jj});
    xlabel('mean ISI, averaged all cells (s)')
    ylabel('probability')
    title(title_label{jj});
    xlim([0 0.1])
end

end