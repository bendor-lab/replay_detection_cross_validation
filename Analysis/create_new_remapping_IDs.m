% CREATE NEW CELL ID FOR REMAPPING STRUCTURE
% Function to convert each session original cell IDs to a new ID in the
% remapping structure.Created to avoid ID repetitions. 
% New ID structure: YXXX - Y = session ID (1 to 5), XXX - orginal cell ID (e.g. 004, 027, 112)
% INPUT:
%      - method: 'wcorr' (default) or 'spearman' 
%      - session: string or cell array - session folder name (from folders_to_process_remapping). If empty, loads general remapping structure
%           with all sessions merged.
%      - original_cell_ID: vector of cell IDs, or cell array of vectors. If empty, loads general remapping structure
%           with all sessions merged.

% we want the option of input: [cell array of folder names and corresponding cell array of cell IDs]
 %                                              [cell array of folder names and vector of cell ids]
 %                                              [string folder name and cell id vector]
 
 function new_ID = create_new_remapping_IDs(method,session,original_cell_ID)
 
 if isempty(original_cell_ID) && isempty(session)
     
     if isempty(method)
         load('X:\BendorLab\Drobo\Neural and Behavioural Data\Rate remapping\Data\rate_remapping_analysis_TRACK_PAIRS_wcorr.mat')
     else
         load(['X:\BendorLab\Drobo\Neural and Behavioural Data\Rate remapping\Data\rate_remapping_analysis_TRACK_PAIRS_' method '.mat']);
     end    
     for s = 1 : length(unique(remapping(1).experiment))
         new_ID = [new_ID; remapping(1).common_good_cells(remapping(1).experiment == s) + (1000 * s)];
     end
     
 elseif ~ (isempty(session) && isempty(original_cell_ID))
    
    load('P:\ground_truth_replay_analysis\Dropbo_data8\folders_to_process_remapping.mat')
    if size(session,1) > 1
        idx_sess = cell2mat(arrayfun(@(x) contains(session,folders{x})'*x, 1:length(folders),'UniformOutput',0)); % index session ID
    else
        idx_sess =  find(arrayfun(@(x) contains(session,folders{x}), 1:length(folders))); % index session ID
    end
    idx_sess(idx_sess == 0) = []; %remove all zeros
    if ischar(session)  % if only one string and a vector of cell IDs convert
        session= {session}; 
    end
    new_ID = [];
    for j=1:length(session)
        if length(strfind(session{j},'\')) > 1 %if session name is the whole path, take last part
            indx = strfind(session{j},'\');
            session{j} = session{j}(indx(end)+1:end);
        end
        % find folder index number in remapping structure
        %idx = find(strcmp(remapping(1).folder, session{j}));
        if length(idx_sess) > 1
            idx = idx_sess(j);
        else
            idx = idx_sess;
        end
        if iscell(original_cell_ID)
            new_ID{j} = original_cell_ID{j} + (1000 * idx);
        elseif length(session)==1 && length(original_cell_ID) >1
            new_ID = original_cell_ID + (1000 * idx);
        else
            new_ID(j) = original_cell_ID(j) + (1000 * idx);
        end
    end
    % for ease of use if only one folder reconvert to array
%     if length(session)==1 && ~iscell(original_cell_ID)
%         new_ID= new_ID{1};
%     end
else 
    error('error in inputs');
end


end