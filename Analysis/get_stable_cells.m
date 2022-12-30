function [experiment,stable_cell_ids,stable_cells_new_ids]= get_stable_cells(folders,method)
% calculates place fields for first and second half of laps for each track
% cells that have peaks > 1Hz in both are kept, rest are discarded
% OUTPUT: indices fo stable cells, coded by experiment, session cell id, and generalised "new id" used in rate_remapping analyses 

load([pwd '\rate_remapping_analysis_TRACK_PAIRS_' method '.mat']);

%initialise
experiment= [];
stable_cell_ids= [];
stable_cells_new_ids= [];

parameters= list_of_parameters;
master_folder= pwd;
for this_folder= 1:length(folders)
    cd([master_folder '\' folders{this_folder}]);
    
    load('lap_times.mat');
    % get first and second lap halves
    midpoint= arrayfun(@(x) round(lap_times(x).total_number_of_laps/2), 1:length(lap_times));
    FIRST_half_ts= arrayfun(@(x) [lap_times(x).start(1) lap_times(x).start(midpoint(x))], 1:length(lap_times),'UniformOutput',0);
    SECOND_half_ts= arrayfun(@(x) [lap_times(x).start(midpoint(x)+1) lap_times(x).end(end)], 1:length(lap_times),'UniformOutput',0);
    
    place_fields_FIRST_half= calculate_place_fields_epochs(parameters.x_bins_width_bayesian,FIRST_half_ts);
    place_fields_SECOND_half= calculate_place_fields_epochs(parameters.x_bins_width_bayesian,SECOND_half_ts);
    
    % get stable cells id for each track
    stable_cells_all= arrayfun(@(x) find(place_fields_FIRST_half.track(x).peak > 1 & place_fields_SECOND_half.track(x).peak > 1), 1:min(length(lap_times),2),'UniformOutput',0);
    % truly stable, common cells need to have peak FR >1 on all 4
    stable_cells_all= intersect(stable_cells_all{1},stable_cells_all{2});
    
    folder_idx= find(strcmp(remapping(2).folder,folders{this_folder}));
    % intersect with cells we are interested in
    cells_tmp= intersect(stable_cells_all,remapping(2).common_good_cells(remapping(2).experiment_all == folder_idx));
    stable_cell_ids= [stable_cell_ids; cells_tmp];
    experiment= [experiment; repmat(folders(this_folder),length(cells_tmp),1)];
    stable_cells_new_ids=  [stable_cells_new_ids; create_new_remapping_IDs('wcorr',folders{this_folder},cells_tmp)];
end
 cd(master_folder);
 
 
if exist('.\Tables\subsets_of_cells.mat')
    load('.\Tables\subsets_of_cells.mat');
    start_row= height(subset_of_cells);
else
    start_row=0;
    subset_of_cells= table(cell(1,1),cell(1,1),cell(1,1),'VariableNames',{'subset','cell_IDs','cdata'});
end

subset_of_cells.subset{start_row+1}= 'stable cells laps';
subset_of_cells.cell_IDs{start_row+1}= stable_cells_new_ids;
save('.\Tables\subsets_of_cells.mat','subset_of_cells');

end