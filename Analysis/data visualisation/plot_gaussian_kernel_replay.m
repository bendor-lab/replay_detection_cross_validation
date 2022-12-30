function plot_gaussian_kernel_replay(folder)
cd folder
load('significant_replay_events_wcorr.mat');

replay_idx= randi(length(significant_replay_events.track(1).spikes));
replay_event= significant_replay_events.track(1).spikes{replay_idx};
replay_dur= significant_replay_events.track(1).event_duration(replay_idx);
replay_time= significant_replay_events.track(1).event_times(replay_idx);
n_sigma= 0.25; 

gauss_ker= gausswin(n_sigma/(1/1000));
gauss_ker= gauss_ker./sum(gauss_ker); 
% we add some time stamps for filtering ease
ts_event= min(replay_time-(replay_dur/2),min(replay_event(:,2)))-n_sigma :1/1000: max(replay_time+(replay_dur/2),max(replay_event(:,2)))+n_sigma;
unique_cells= unique(replay_event(:,1));

figure('Color','w'); 
ax1= subplot(2,1,1); hold on;
arrayfun(@(x) plot(replay_event(replay_event(:,1)==unique_cells(x),2),x*ones(size(replay_event(replay_event(:,1)==unique_cells(x),1))),'*'),...
    1:length(unique_cells));
ax2= subplot(2,1,2); hold on;
for this_cell=1:length(unique_cells)
    cell_spikes= replay_event(replay_event(:,1)==unique_cells(this_cell),2);
    example_spike_train= histcounts(cell_spikes,ts_event);
    filt_sig = filter(gauss_ker,1,example_spike_train);
    max_FR(this_cell)= max(filt_sig);
    plot(ts_event(1:end-1)+0.5*mean(diff(ts_event)),filt_sig);
end
title(['filter length ' num2str(n_sigma*1000) ' ms, sigma= ' num2str((n_sigma*1000-1)/(2*2.5)) ' ms'])
linkaxes([ax1 ax2],'x')


% SCALING
scale_factor= ones(length(unique_cells),1);
scale_factor(1)= max_FR(2)./max_FR(1);

figure('Color','w'); 
ax1= subplot(2,1,1); hold on;
arrayfun(@(x) plot(replay_event(replay_event(:,1)==unique_cells(x),2),x*ones(size(replay_event(replay_event(:,1)==unique_cells(x),1))),'*'),...
    1:length(unique_cells));
ax2= subplot(2,1,2); hold on;
for this_cell=1:length(unique_cells)
    cell_spikes= replay_event(replay_event(:,1)==unique_cells(this_cell),2);
    example_spike_train= histcounts(cell_spikes,ts_event)*scale_factor(this_cell);
    filt_sig = filter(gauss_ker,1,example_spike_train);
    max_FR(this_cell)= max(filt_sig);
    plot(ts_event(1:end-1)+0.5*mean(diff(ts_event)),filt_sig);
end
title(['filter length ' num2str(n_sigma*1000) ' ms, sigma= ' num2str((n_sigma*1000-1)/(2*2.5)) ' ms'])
linkaxes([ax1 ax2],'x')



end