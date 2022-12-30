function prepare_data_for_regression(varargin)
% INPUTS: name value pairs
%               'option': 'wcorr','spearman' or 'linear' (in the future)
%               'epoch': any epoch in rate_remapping_TRACK_PAIRS_option.mat
% OUTPUT:
% saves data_for_regression table of cell ID - event - track replayed- num spikes etc...
% under .\Tables\replay_spiking_data_regression.mat

p= inputParser;
addParameter(p,'option','wcorr',@ischar);
addParameter(p,'epoch','POST',@ismatrix);
parse(p,varargin{:});
option= p.Results.option;

load(['rate_remapping_analysis_TRACK_PAIRS_' option '.mat']);
row_idx= find(strcmp([remapping_raw.epoch],p.Results.epoch));

data_for_regression= table;
observations= 1;
% new_cell_id=1;
warning('off');

% we purposely create a table in the long format, for R
for this_session= 1:length(unique(remapping_raw(row_idx).experiment))
    cell_ids= remapping_raw(row_idx).new_ID(remapping_raw(row_idx).experiment==this_session);
    
    for this_cell= 1:size(remapping_raw(row_idx).replay_spikes1{this_session},1)
        % start with track 1 replay events
        non_zero_evts1= remapping_raw(row_idx).replay_spikes1{this_session}(this_cell,:)~=0;
        non_zero_evts2= remapping_raw(row_idx).replay_spikes2{this_session}(this_cell,:)~=0;
        for this_event= 1:size(remapping_raw(row_idx).replay_spikes1{this_session},2)
            if remapping_raw(row_idx).replay_spikes1{this_session}(this_cell,this_event)~=0
                data_for_regression.cell_ID(observations)= cell_ids(this_cell);
                data_for_regression.event(observations)= this_event; % probably not useful
                data_for_regression.track_replayed(observations)= 1;
%                 data_for_regression.num_spikes(observations)= remapping_raw(row_idx).replay_spikes1{this_session}(this_cell,this_event);
                data_for_regression.rate(observations)= remapping_raw(row_idx).max_instantaneous_FR_1{this_session}(this_cell,this_event);
                % need to change to get zscore across both types events not
                % one track only
                data_for_regression.zscored_rate(observations)= (remapping_raw(row_idx).max_instantaneous_FR_1{this_session}(this_cell,this_event) ...
                                                                                                            - mean([remapping_raw(row_idx).max_instantaneous_FR_1{this_session}(this_cell,non_zero_evts1)...
                                                                                                            remapping_raw(row_idx).max_instantaneous_FR_2{this_session}(this_cell,non_zero_evts2)]))...
                                                                                                            ./std([remapping_raw(row_idx).max_instantaneous_FR_1{this_session}(this_cell,non_zero_evts1)...
                                                                                                            remapping_raw(row_idx).max_instantaneous_FR_2{this_session}(this_cell,non_zero_evts2)]);
                observations= observations+1;
            end
        end       
        % then track 2 replay events
        non_zero_evts1= remapping_raw(row_idx).replay_spikes1{this_session}(this_cell,:)~=0;
        non_zero_evts2= remapping_raw(row_idx).replay_spikes2{this_session}(this_cell,:)~=0;
        for this_event= 1:size(remapping_raw(row_idx).replay_spikes2{this_session},2)
            if remapping_raw(row_idx).replay_spikes2{this_session}(this_cell,this_event)~=0
                data_for_regression.cell_ID(observations)= cell_ids(this_cell);
                data_for_regression.event(observations)= this_event; % probably not useful
                data_for_regression.track_replayed(observations)= 2;
%                 data_for_regression.num_spikes(observations)= remapping_raw(row_idx).replay_spikes2{this_session}(this_cell,this_event);
                data_for_regression.rate(observations)= remapping_raw(row_idx).max_instantaneous_FR_2{this_session}(this_cell,this_event);
                data_for_regression.zscored_rate(observations)= (remapping_raw(row_idx).max_instantaneous_FR_2{this_session}(this_cell,this_event) ...
                                                                                                       - mean([remapping_raw(row_idx).max_instantaneous_FR_1{this_session}(this_cell,non_zero_evts1)...
                                                                                                        remapping_raw(row_idx).max_instantaneous_FR_2{this_session}(this_cell,non_zero_evts2)]))...
                                                                                                        ./std([remapping_raw(row_idx).max_instantaneous_FR_1{this_session}(this_cell,non_zero_evts1)...
                                                                                                        remapping_raw(row_idx).max_instantaneous_FR_2{this_session}(this_cell,non_zero_evts2)]);
                observations= observations+1;
            end
        end       
%         new_cell_id= new_cell_id+1; % this way each cell has its unique identifier (not same cells across sessions)
    end
end

save('./Tables/replay_spiking_data_regression.mat','data_for_regression');
% writetable(data_for_regression,'./Tables/replay_spiking_data_regression.txt')

%% analysis bit 
warning('off');
unique_cells= unique(data_for_regression.cell_ID);
data_for_regression.track_replayed= data_for_regression.track_replayed>1;
spikes_edges= [0.5:1:10.5]; spike_ctrs= [1:1:10];
rate_edges= [0:5:70]; rate_ctrs= [0.5:5:69.5];
figure('Color','w');
pval=[]; weight=[]; num_events=[];
for this_cell=1:length(unique_cells)
    m= fitglm(data_for_regression(data_for_regression.cell_ID== unique_cells(this_cell),:),...
       'track_replayed ~ zscored_rate','Distribution','binomial');
   weight(this_cell)= m.Coefficients.Estimate(2);
   pval(this_cell)= m.Coefficients.pValue(2);
   num_events(this_cell)= length(find(data_for_regression.cell_ID== unique_cells(this_cell)));
   
   subplot(ceil(sqrt(length(unique_cells))),ceil(sqrt(length(unique_cells))),this_cell)
   if pval(this_cell) <= 0.05
       set(gca,'Color','c');
   end
    n1= histcounts(data_for_regression.rate(data_for_regression.cell_ID==unique_cells(this_cell) &...
                                data_for_regression.track_replayed==0),rate_edges,'Normalization','probability');
    n2= histcounts(data_for_regression.rate(data_for_regression.cell_ID==unique_cells(this_cell) &...
                                data_for_regression.track_replayed==1),rate_edges,'Normalization','probability');
    hold on; plot(rate_ctrs,n1,'k'); plot(rate_ctrs,n2,'r');
    title(['Cell ' num2str(unique_cells(this_cell))])   
end

cdata= repmat([0.7 0.7 0.7],length(unique_cells),1);
predictive_cells= unique_cells(pval <= 0.05);
cdata(pval<=0.05,1)=1;cdata(pval<=0.05,2:3)=0;

fig2= figure('Color','w');
ax1= subplot(1,3,1); hold on;
arrayfun(@(x) plot(sort([0 weight(x)]),[x x],'Color',cdata(x,:),'LineWidth',1.5), 1:length(weight));
% xlim([min(weight)-1 max(weight)]); 
ylim([0 length(weight)+1]);
xlim([-4 4]); ylabel(ax1,'cell ID'); xlabel(ax1,'regression weight')
title(ax1,['percentage of predictive cells: ' num2str(100*length(predictive_cells)/length(unique_cells)) '%'])
ax2= subplot(1,3,2); hold on;
arrayfun(@(x) plot(sort([0 pval(x)]),[x x],'Color',cdata(x,:),'LineWidth',1.5), 1:length(pval));
xlim([-0.005 0.2]); ylim([0 length(pval)+1]);
plot([0.05 0.05],ylim,'k--','LineWidth',1.5)
 xlabel(ax2,'p-value')
% title([{'zscored rate'}; {'+ logistic regression each cell'}])
subplot(1,3,3); hold on;
arrayfun(@(x) plot(sort([0 num_events(x)]),[x x],'Color',cdata(x,:),'LineWidth',1.5), 1:length(num_events));
ylim([0 length(pval)+1]);
xlabel([{'number of replay events';'cell participates in'}])
run_format_settings(fig2)

if exist('.\Tables\subsets_of_cells.mat')
    load('.\Tables\subsets_of_cells.mat');
    start_row= height(subset_of_cells);
else
    start_row=0;
    subset_of_cells= table(cell(1,1),cell(1,1),cell(1,1),'VariableNames',{'subset','cell_IDs','cdata'});
end

subset_of_cells.subset{start_row+1}= 'predictive cells';
subset_of_cells.cell_IDs{start_row+1}= predictive_cells;

save('.\Tables\subsets_of_cells.mat','subset_of_cells');
warning('on');

end