% WITHIN FIELD ANALYSIS
% Code to extract area under the curve (AUC) for each place field, and to
% detect place cells with multiple fields
% INPUT:
    % method: 'wcorr'(default) or 'spearman'. 

function within_field_analysis(method)

load('X:\BendorLab\Drobo\Neural and Behavioural Data\Rate remapping\Data\folders_to_process_remapping.mat')
if strcmp(method,'spearman')
    load('X:\BendorLab\Drobo\Neural and Behavioural Data\Rate remapping\Data\rate_remapping_analysis_TRACK_PAIRS_spearman.mat')
else
   load('X:\BendorLab\Drobo\Neural and Behavioural Data\Rate remapping\Data\rate_remapping_analysis_TRACK_PAIRS_wcorr.mat')
end

for fol = 1 : length(folders)
    
    load([folders{fol} '\extracted_place_fields_BAYESIAN.mat'])
    ct = 1;
    plfield = [];
    
    % Get common cells from both tracks
    common_good_cells= intersect(place_fields_BAYESIAN.track(1).good_cells, place_fields_BAYESIAN.track(2).good_cells);

    for track = 1 : 2
        
        % For each good cell in the track
        for pc = 1 : length(common_good_cells)
            pc_bounds = [];
            final_pc_bounds = [];
            cell_id = common_good_cells(pc);
            plfield(ct).track = track;
            plfield(ct).cell_id = cell_id;
            
            %figure;plot(place_fields.track(track).raw{cell_id})
            
            % Define place field as region where firing rate goes above 20 of max
            [max_firing, max_idx] = max(place_fields_BAYESIAN.track(track).raw{cell_id});
            baseline_firing = min(place_fields_BAYESIAN.track(track).raw{cell_id});
            firing_thresh = baseline_firing +(0.2 * (max_firing - baseline_firing));
            firing_idx = place_fields_BAYESIAN.track(track).raw{cell_id} > firing_thresh;
            edge_detect_idx = diff(firing_idx);
            
            %hold on
            %plot([min(xlim) max(xlim)],[firing_thresh firing_thresh],'--')
            
            % Check if place field is only one bin in the edge
            indices = edge_detect_idx(find(edge_detect_idx ~= 0));
            edge_idx = find(edge_detect_idx ~= 0);
            if  length(indices) == 1 & edge_idx == 1
                final_pc_bounds = [1 2];
            elseif length(indices) == 1 & edge_idx == length(edge_detect_idx)
                final_pc_bounds = [length(edge_detect_idx)-1 length(edge_detect_idx)]; %19 & 20 in 200cm track with 10cm bins
            else
                % Check if edges are around max peak
                if ~any(find(edge_detect_idx == 1) < max_idx) |  ~any(find(edge_detect_idx == -1) > max_idx)
                    indices = edge_detect_idx(edge_detect_idx ~= 0);
                    if indices(1) == -1
                        edge_detect_idx(1) = 1; % means that FR starts above 0
                    end
                    if indices(end) == 1
                        edge_detect_idx(end) = -1; %means that FR ends above 0
                    end
                end
                
                %if odd, means that there's a indx missing
                if mod(length(find(edge_detect_idx ~= 0)),2) ~= 0
                    diff_indices = diff(edge_detect_idx(edge_detect_idx ~= 0));
                    if diff_indices(end) ==  0 & edge_detect_idx(end) == -1 % if indices end as -1 -1 and it's in the edge
                        pc_bounds = [length(edge_detect_idx)-1 length(edge_detect_idx)];
                    elseif diff_indices(1) ==  0 & edge_detect_idx(end) == 1 % if indices start as 1 1 and it's in the edge
                        pc_bounds = [1 2];
                    else
                        if indices(1) == -1
                            edge_detect_idx(1) = 1; % means that FR starts above 0
                        end
                        if indices(end) == 1
                            edge_detect_idx(end) = -1; %means that FR ends above 0
                        end
                    end
                end
                % Find place cell boundaries
                if isempty(pc_bounds)
                    pc_bounds(:,1) = find(edge_detect_idx == 1);
                    % Check again if there are different number of start and end boundaries
                    if length(find(edge_detect_idx == -1)+1) ~= length(pc_bounds(:,1))
                        idx1 = find(diff(edge_detect_idx(edge_detect_idx ~=0)) == 0);
                        idcs = find(edge_detect_idx ~=0);
                        edge_detect_idx(idcs(idx1)) = 0;
                        pc_bounds = [];
                        pc_bounds(:,1) = find(edge_detect_idx == 1);
                        pc_bounds(:,2)= find(edge_detect_idx == -1)+1;
                    else
                        pc_bounds(:,2)= find(edge_detect_idx == -1)+1; % to account for diff
                    end
                    % Check that end boundary is not before start boundary
                    if pc_bounds(:,1) > pc_bounds(:,2)
                        pc_bounds = [];
                        edge_detect_idx(1) = 1; % means that FR starts above 0
                        edge_detect_idx(end) = -1; %means that FR ends above 0
                        if length(find(edge_detect_idx == -1)+1) ~= length(find(edge_detect_idx == 1)) %if there are different number of start and end boundaries
                            idx1 = find(diff(edge_detect_idx(edge_detect_idx ~=0)) == 0);
                            idcs = find(edge_detect_idx ~=0);
                            edge_detect_idx(idcs(idx1)) = 0;
                            pc_bounds = [];
                            pc_bounds(:,1) = find(edge_detect_idx == 1);
                            pc_bounds(:,2)= find(edge_detect_idx == -1)+1;
                        else
                            pc_bounds(:,1) = find(edge_detect_idx == 1);
                            pc_bounds(:,2)= find(edge_detect_idx == -1)+1;
                        end
                    end
                    % Find threshold crossing that includes max (in case threshold is crossed multiple times)
                    pf_id = find(pc_bounds(:,2) >= max_idx,1,'first'); % choose the boundaries of the peak with the max FR
                    final_pc_bounds = pc_bounds(pf_id,:);
                else
                    final_pc_bounds = pc_bounds;
                end
            end
            
            % Calculate and save area under the curve for that place field
            plfield(ct).AUC = trapz(place_fields_BAYESIAN.track(track).raw{cell_id}(final_pc_bounds(1):final_pc_bounds(2))- baseline_firing);
            
            %hold on
            %area(final_pc_bounds(1):final_pc_bounds(2),[place_fields.track(track).raw{cell_id}(final_pc_bounds(1):final_pc_bounds(2)) - baseline_firing])
            
            % Check if there are extra peaks
            if size(pc_bounds,1) > 1
                other_pc_bounds = pc_bounds;
                other_pc_bounds(pc_bounds == final_pc_bounds) =[];
                other_pc_bounds = reshape(other_pc_bounds,[length(other_pc_bounds)/2,2]); %reshape
                if (other_pc_bounds(:,2) - other_pc_bounds(:,1)) <= 0 % exceptions with very noisy cells
                    ct = ct+1;
                    continue
                end
                % Finds extra peaks FR and the corresponding indices
                for k = 1 : size(other_pc_bounds,1)
                    extra_peaks(k) = max(place_fields_BAYESIAN.track(track).raw{cell_id}(other_pc_bounds(k,1):other_pc_bounds(k,2)));
                    idcs = find(place_fields_BAYESIAN.track(track).raw{cell_id} == extra_peaks(k));
                    if length(idcs) > 1
                        [~,d] = min(abs(other_pc_bounds(k,:) - idcs));
                        extra_peaks_idx(k) = idcs(d);
                    else
                        extra_peaks_idx(k) = idcs;
                    end
                end
                % check that the extra peaks do not overlap with the MAX
                % peak, that cover more than 2 bins, and the peak distance
                % is more than 4 bins (40cm)
                if (other_pc_bounds > final_pc_bounds |  other_pc_bounds < final_pc_bounds) & any(other_pc_bounds(:,2)- other_pc_bounds(:,1) > 1)...
                        & any(abs(max_idx - extra_peaks_idx) > 4)
                    plfield(ct).multi_plfields = 1;
                end
                clear extra_peaks extra_peaks_idx
            end
            ct = ct+1;
        end
    end
    
    session(fol).plfield = plfield;
    %close all
end

% SAVE IN REMAPPING STRUCTURE
for per = 1 : size(remapping,1) % for each epoch 
    ct = 1;
    for f = 1 : length(folders) % for each session
        AUCS = session(f).plfield;
        idx = remapping(per).experiment == f;
        cells_id = remapping(per).ID_active_cells_during_replay(idx);
        for ce = 1 : length(cells_id) %for each common cell active in replay
            t1 = AUCS([AUCS(:).track] == 1);
            t2 = AUCS([AUCS(:).track] == 2);
            T1_AUC = t1([t1(:).cell_id] == cells_id(ce)).AUC;
            T2_AUC = t2([t2(:).cell_id] == cells_id(ce)).AUC;
            remapping_raw(per).plfield_AUC_1(ct) = T1_AUC;
            remapping_raw(per).plfield_AUC_2(ct) = T2_AUC;
            remapping(per).plfield_AUC_diff(ct) = T1_AUC-T2_AUC;
            % Check if the cell has multiple place fields in any of the tracks
            if t1([t1(:).cell_id] == cells_id(ce)).multi_plfields == 1 | t2([t2(:).cell_id] == cells_id(ce)).multi_plfields == 1
                remapping(per).multi_plfield_cellID(ct) = cells_id(ce) + (1000*f); %  save as new cell ID
                if t1([t1(:).cell_id] == cells_id(ce)).multi_plfields == 1
                    remapping_raw(per).T1_multi_plfield_cellID(ct) = cells_id(ce) + (1000*f);
                end
                if  t2([t2(:).cell_id] == cells_id(ce)).multi_plfields == 1
                    remapping_raw(per).T2_multi_plfield_cellID(ct) = cells_id(ce) + (1000*f);
                end
            end
            ct = ct + 1;
            
        end
    end
    remapping(per).multi_plfield_cellID(remapping(per).multi_plfield_cellID == 0) = [];
    remapping_raw(per).T1_multi_plfield_cellID(remapping_raw(per).T1_multi_plfield_cellID == 0) = [];
    remapping_raw(per).T2_multi_plfield_cellID(remapping_raw(per).T2_multi_plfield_cellID == 0) = [];
end

% save in subset_of_cells a vector with cell IDs of cells with a single place field
load('.\Tables\subsets_of_cells.mat');
temp = [];
for per = 1 : size(remapping,1) % for each epoch
    temp = [temp remapping(per).new_ID(~ismember(remapping(per).new_ID,remapping(per).multi_plfield_cellID))];
end
if any(strcmp(subset_of_cells.subset,'single place field cells')) %if already exists
    row_id = find(strcmp(subset_of_cells.subset,'single place field cells'));
    subset_of_cells.cell_IDs{row_id} = unique(temp);
else
    new_row = size(subset_of_cells,1)+1;
    subset_of_cells.subset{new_row} = 'single place field cells';
    subset_of_cells.cell_IDs{new_row} = unique(temp);
end
if any(strcmp(subset_of_cells.subset,'multiple place field cells')) %if already exists
    row_id = find(strcmp(subset_of_cells.subset,'multiple place field cells'));
    subset_of_cells.cell_IDs{row_id} = unique([remapping(:).multi_plfield_cellID]);
else
    new_row = size(subset_of_cells,1)+1;
    subset_of_cells.subset{new_row} = 'multiple place field cells';
    subset_of_cells.cell_IDs{new_row} = unique([remapping(:).multi_plfield_cellID]);
end

cd('X:\BendorLab\Drobo\Neural and Behavioural Data\Rate remapping\Data')        
save multiple_place_fields_all_sessions_wcorr session
save('.\Tables\subsets_of_cells.mat','subset_of_cells')
switch method
    case 'wcorr'
        save rate_remapping_analysis_TRACK_PAIRS_wcorr remapping remapping_raw
        %save('.\Tables\subsets_of_cells_wcorr.mat','subset_of_cells')
    case 'spearman'
        save rate_remapping_analysis_TRACK_PAIRS_spearman remapping remapping_raw
        %save('.\Tables\subsets_of_cells_spearman.mat','subset_of_cells')
    otherwise
        save rate_remapping_analysis_TRACK_PAIRS remapping remapping_raw
end

end
                