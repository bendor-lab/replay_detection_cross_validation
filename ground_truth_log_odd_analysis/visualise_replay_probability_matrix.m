
function [] = visualise_replay_probability_matrix(log_odd,folders,folder_number)

index_pval = intersect(find(log_odd.experiment==folder_number),...
    find(log_odd.pvalue >= log10(0.001) & log_odd.pvalue <= log10(0.005)));
index_log_odd = intersect(index_pval,find(...
    log_odd.normal_zscored.original>= -1 & ...
    log_odd.normal_zscored.original <= 3.5));

% index_pval = intersect(log_odd.index(log_odd.experiment==folder_number),...
%     log_odd.index(log_odd.pvalue >= log10(0.001) & log_odd.pvalue <= log10(0.005)));
% index_log_odd = intersect(index_pval,log_odd.index(...
%     log_odd.normal_zscored.original<= 0.2 & ...
%     log_odd.normal_zscored.original >= -0.2));


event_index = index_log_odd;
index = log_odd.index(event_index);

event_index(20)

%
% index = intersect(index_log_odd,log_odd.index(log_odd.behavioural_state == 0));

% index = intersect(log_odd.index(log_odd.experiment==folder_number),log_odd.index(log_odd.pvalue<= log10(0.01)));
% index = [799 1008 1011 1067 1130 1154 1244 1268 1300 1636];

cd(folders{folder_number})
load estimated_position_20ms
nfig = 1;
limits = [];

for event = ceil(1:2:length(index));
    
    figure(nfig)
    subplot(2,2,1)
    imagesc(estimated_position_20ms(1).replay_events(index(event)).replay(1:20,:));
    colormap(flipud(bone))
    colorbar
    title({sprintf('Track %i event %i \nz-scored log odd = %.2f & p-value =%.2f ',...
        log_odd.track(event_index(event)),index(event),log_odd.normal_zscored.original(event_index(event)),log_odd.pvalue(event_index(event)))})
    max_point = max(max(estimated_position_20ms(1).replay_events(index(event)).replay));
    min_point = min(min(estimated_position_20ms(1).replay_events(index(event)).replay));
    xlabel('Time bin')
    ylabel('Track 1 Position bin')
    limits = [min_point, max_point];
    caxis(limits)
    
    subplot(2,2,2)
    imagesc(estimated_position_20ms(1).replay_events(index(event)).replay(21:40,:));
    colormap(flipud(bone))
    colorbar
    title({sprintf('Track %i event %i \nz-scored log odd = %.2f & p-value =%.2f ',...
        log_odd.track(event_index(event)),index(event),log_odd.normal_zscored.original(event_index(event)),log_odd.pvalue(event_index(event)))})
    max_point = max(max(estimated_position_20ms(1).replay_events(index(event)).replay));
    min_point = min(min(estimated_position_20ms(1).replay_events(index(event)).replay));
    xlabel('Time bin')
    ylabel('Track 1 Position bin')
    limits = [min_point, max_point];
    caxis(limits)
    
    subplot(2,2,3)
    imagesc(estimated_position_20ms(1).replay_events(index(event+1)).replay(1:20,:));
    colormap(flipud(bone))
    colorbar
    title({sprintf('Track %i event %i \nz-scored log odd = %.2f & p-value =%.2f ',...
        log_odd.track(event_index(event)),index(event),log_odd.normal_zscored.original(event_index(event)),log_odd.pvalue(event_index(event)))})
    max_point = max(max(estimated_position_20ms(1).replay_events(index(event+1)).replay));
    min_point = min(min(estimated_position_20ms(1).replay_events(index(event+1)).replay));
    xlabel('Time bin')
    ylabel('Track 1 Position bin')
    limits = [min_point, max_point];
    caxis(limits)
    
    subplot(2,2,4)
    imagesc(estimated_position_20ms(1).replay_events(index(event+1)).replay(21:40,:));
    colormap(flipud(bone))
    colorbar
    title({sprintf('Track %i event %i \nz-scored log odd = %.2f & p-value =%.2f ',...
        log_odd.track(event_index(event)),index(event),log_odd.normal_zscored.original(event_index(event)),log_odd.pvalue(event_index(event)))})
    max_point = max(max(estimated_position_20ms(1).replay_events(index(event+1)).replay));
    min_point = min(min(estimated_position_20ms(1).replay_events(index(event+1)).replay));
    xlabel('Time bin')
    ylabel('Track 1 Position bin')
    limits = [min_point, max_point];
    caxis(limits)
    
    
    nfig = nfig +1
    sgtitle('Posterior Probability after decoding')
end
cd ..