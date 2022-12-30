function [remapping, remapping_raw] = rate_remapping_TRACK_PAIRS(folders,varargin)
% Input:
    % folders is [] for individual session, or cell array of folders
    % varargin{1} is for methods, and can be 'spearman' or 'wcorr'
    % varargin{2} as save toggle - can be 1 (default) or 0
% Output: 
    % remapping: difference between variables
    % remapping_raw: variables for each track

if isempty(folders)
%     load('X:\BendorLab\Drobo\Neural and Behavioural Data\Rate remapping\Data\folders_to_process_remapping.mat')
    a=pwd;
    folders={a};
end
master_folder= pwd;

if ~isempty(varargin)
    method= varargin{1};
    if length(varargin)>1
        save_option = varargin{2};
    else
        save_option= 1;
    end
    disp(['running for... ' varargin{1}])
else
    method= 'wcorr';
    save_option = 1;
end


% Allocate variables in two structures
track_pair = 1; %not removing but would only be used if we ever have more than two tracks
for epoch = 1 : 7 % PRE & POST, sleep and awake   
    remapping(epoch,track_pair).experiment=[];
    remapping(epoch,track_pair).experiment_all= [];
    remapping(epoch,track_pair).folder= [];
    remapping(epoch,track_pair).epoch= [];
    remapping(epoch,track_pair).common_good_cells=[];
    remapping(epoch,track_pair).ID_active_cells_during_replay=[];
    remapping(epoch,track_pair).new_ID =[];
    remapping(epoch,track_pair).PRE_to_POST_active_cells=[];
    remapping(epoch,track_pair).fraction_of_active_cells_during_replay=[];
    remapping(epoch,track_pair).place_field_diff=[];
    remapping(epoch,track_pair).log_place_field_diff= [];
    remapping(epoch,track_pair).place_field_centre_diff=[];
    remapping(epoch,track_pair).replay_spike_diff=[];
    remapping(epoch,track_pair).replay_spike_diff_nonZero=[];
    remapping(epoch,track_pair).replay_spike_median_diff_nonZero=[];
    remapping(epoch,track_pair).replay_rate_diff=[];
    remapping(epoch,track_pair).log_replay_rate_diff=[];
    remapping(epoch,track_pair).median_replay_rate_diff=[];
    remapping(epoch,track_pair).track_proportion_events_diff=[];
    remapping(epoch,track_pair).track_mean_rate_diff= [];
    remapping(epoch,track_pair).mean_max_FR_replay_diff= [];
    
    remapping_raw(epoch,track_pair).experiment=[];
    remapping_raw(epoch,track_pair).folder= [];
    remapping_raw(epoch,track_pair).epoch= [];
    remapping_raw(epoch,track_pair).common_good_cells=[];
    remapping_raw(epoch,track_pair).ID_active_cells_during_replay=[];
    remapping_raw(epoch,track_pair).new_ID =[];
    remapping_raw(epoch,track_pair).PRE_to_POST_active_cells=[];
    remapping_raw(epoch,track_pair).raw_peak_BAYESIAN_plfield_1=[];
    remapping_raw(epoch,track_pair).raw_peak_BAYESIAN_plfield_2=[];
    remapping_raw(epoch,track_pair).track1_mean_replay_spikes=[];
    remapping_raw(epoch,track_pair).track2_mean_replay_spikes=[];
    remapping_raw(epoch,track_pair).track1_mean_replay_spikes_nonZero=[];
    remapping_raw(epoch,track_pair).track2_mean_replay_spikes_nonZero=[];
    remapping_raw(epoch,track_pair).track1_median_replay_spikes_nonZero=[];
    remapping_raw(epoch,track_pair).track2_median_replay_spikes_nonZero=[];
    remapping_raw(epoch,track_pair).track1_median_replay_rate=[];
    remapping_raw(epoch,track_pair).track2_median_replay_rate=[];
    remapping_raw(epoch,track_pair).replay_spikes1=[];
    remapping_raw(epoch,track_pair).replay_spikes2=[];
    remapping_raw(epoch,track_pair).replay_rate1=[];
    remapping_raw(epoch,track_pair).replay_rate2=[];
    remapping_raw(epoch,track_pair).proportion_events_cell_active_track1= [];
    remapping_raw(epoch,track_pair).proportion_events_cell_active_track2= [];
    remapping_raw(epoch,track_pair).mean_rate_T1= [];
    remapping_raw(epoch,track_pair).mean_rate_T2= [];
    remapping_raw(epoch,track_pair).proportion_events_T1= [];
    remapping_raw(epoch,track_pair).proportion_events_T2= [];
    remapping_raw(epoch,track_pair).max_instantaneous_FR_1= [];
    remapping_raw(epoch,track_pair).max_instantaneous_FR_2= [];
    remapping_raw(epoch,track_pair).track1_mean_replay_inst_FR_nonZero = [];
    remapping_raw(epoch,track_pair).track2_mean_replay_inst_FR_nonZero = [];
    
%       remapping(epoch,track_pair).peak_FR_T1=[];
%       remapping(epoch,track_pair).peak_FR_T2=[];    
end
cd(master_folder);

for i = 1 : length(folders)
    cd(folders{i});
    sort_replay_events([],'wcorr')
    number_of_significant_replays(0.05,3, 'wcorr', [])
    
    data = compare_replay_across_tracks(method);
    for this_epoch=1:7
        % Get new cell IDs
        if ~isempty(data(this_epoch).ID_active_cells_during_replay) %if there's PRE session
            data(this_epoch).new_ID = create_new_remapping_IDs(method,folders{i},data(this_epoch).ID_active_cells_during_replay);
        end
    end
%     data(2).new_ID = create_new_remapping_IDs(method,folders{i},data(2).ID_active_cells_during_replay);
    % Find common cells that are active from PRE to POST and save
    [data(1).PRE_to_POST_active_cells, data(2).PRE_to_POST_active_cells ] = deal(intersect(data(2).new_ID,data(1).new_ID));
    [data(3).PRE_to_POST_active_cells, data(4).PRE_to_POST_active_cells ] = deal(intersect(data(3).new_ID,data(4).new_ID));
    [data(5).PRE_to_POST_active_cells, data(6).PRE_to_POST_active_cells ] = deal(intersect(data(5).new_ID,data(6).new_ID));

    for epoch = 1 : size(data,1)
        for track_pair = 1 : size(data,2)
            remapping(epoch,track_pair).folder = [remapping(epoch,track_pair).folder folders(i)];
            remapping_raw(epoch,track_pair).folder = [remapping_raw(epoch,track_pair).folder folders(i)];
            if epoch==1
                remapping(epoch,track_pair).epoch= {'sleepPRE'};
                remapping_raw(epoch,track_pair).epoch= {'sleepPRE'};                
            elseif epoch ==2
                remapping(epoch,track_pair).epoch= {'sleepPOST'};
                remapping_raw(epoch,track_pair).epoch= {'sleepPOST'};
            elseif epoch ==3
                remapping(epoch,track_pair).epoch= {'awakePRE'};
                remapping_raw(epoch,track_pair).epoch= {'awakePRE'};
            elseif epoch ==4
                remapping(epoch,track_pair).epoch= {'awakePOST'};
                remapping_raw(epoch,track_pair).epoch= {'awakePOST'};
            elseif epoch==5
                remapping(epoch,track_pair).epoch= {'PRE'};
                remapping_raw(epoch,track_pair).epoch= {'PRE'};
            elseif epoch==6
                remapping(epoch,track_pair).epoch= {'POST'};
                remapping_raw(epoch,track_pair).epoch= {'POST'};
            elseif epoch==7
                remapping(epoch,track_pair).epoch= {'RUN'};
                remapping_raw(epoch,track_pair).epoch= {'RUN'};
            end
            remapping(epoch,track_pair).experiment = [remapping(epoch,track_pair).experiment; (i*ones(size(data(epoch,track_pair).place_field_diff)))];
            remapping_raw(epoch,track_pair).experiment = [remapping_raw(epoch,track_pair).experiment; (i*ones(size(data(epoch,track_pair).place_field_diff)))];
            remapping(epoch,track_pair).experiment_all = [remapping(epoch,track_pair).experiment_all; (i*ones(size(data(epoch,track_pair).good_cells')))];
            remapping(epoch,track_pair).common_good_cells = [remapping(epoch,track_pair).common_good_cells; data(epoch,track_pair).good_cells'];
            remapping_raw(epoch,track_pair).common_good_cells = [remapping_raw(epoch,track_pair).common_good_cells; data(epoch,track_pair).good_cells'];
            remapping(epoch,track_pair).ID_active_cells_during_replay = [remapping(epoch,track_pair).ID_active_cells_during_replay; data(epoch,track_pair).ID_active_cells_during_replay'];
            remapping_raw(epoch,track_pair).ID_active_cells_during_replay = [remapping_raw(epoch,track_pair).ID_active_cells_during_replay; data(epoch,track_pair).ID_active_cells_during_replay'];
            remapping(epoch,track_pair).new_ID = [remapping(epoch,track_pair).new_ID data(epoch,track_pair).new_ID];
            remapping_raw(epoch,track_pair).new_ID = [remapping_raw(epoch,track_pair).new_ID data(epoch,track_pair).new_ID];
            remapping(epoch,track_pair).fraction_of_active_cells_during_replay = [remapping(epoch,track_pair).fraction_of_active_cells_during_replay; data(epoch,track_pair).fraction_of_active_cells_during_replay];
            remapping(epoch,track_pair).PRE_to_POST_active_cells = [remapping(epoch,track_pair).PRE_to_POST_active_cells; data(epoch,track_pair).PRE_to_POST_active_cells'];
            remapping_raw(epoch,track_pair).PRE_to_POST_active_cells = [remapping_raw(epoch,track_pair).PRE_to_POST_active_cells; data(epoch,track_pair).PRE_to_POST_active_cells'];

            % Remapping- diff tracks structure
            remapping(epoch,track_pair).place_field_diff =[remapping(epoch,track_pair).place_field_diff; data(epoch,track_pair).place_field_diff];
            remapping(epoch,track_pair).log_place_field_diff = [remapping(epoch,track_pair).log_place_field_diff; data(epoch,track_pair).place_field_diff];
            remapping(epoch,track_pair).place_field_centre_diff = [remapping(epoch,track_pair).place_field_centre_diff; data(epoch,track_pair).place_field_centre_diff];
            remapping(epoch,track_pair).replay_spike_diff = [remapping(epoch,track_pair).replay_spike_diff; data(epoch,track_pair).replay_spike_diff];
            remapping(epoch,track_pair).replay_spike_diff_nonZero = [remapping(epoch,track_pair).replay_spike_diff_nonZero; data(epoch,track_pair).replay_spike_diff_nonZero];
            remapping(epoch,track_pair).replay_spike_median_diff_nonZero = [remapping(epoch,track_pair).replay_spike_median_diff_nonZero; data(epoch,track_pair).replay_spike_median_diff_nonZero];
            remapping(epoch,track_pair).median_replay_rate_diff = [remapping(epoch,track_pair).median_replay_rate_diff; data(epoch,track_pair).median_replay_rate_diff];
            remapping(epoch,track_pair).replay_rate_diff = [remapping(epoch,track_pair).replay_rate_diff; data(epoch,track_pair).replay_rate_diff];
            remapping(epoch,track_pair).log_replay_rate_diff = [remapping(epoch,track_pair).log_replay_rate_diff; data(epoch,track_pair).log_replay_rate_diff];
            remapping(epoch,track_pair).track_proportion_events_diff = [remapping(epoch,track_pair).track_proportion_events_diff data(epoch,track_pair).proportion_events_diff];
            remapping(epoch,track_pair).track_mean_rate_diff = [remapping(epoch,track_pair).track_mean_rate_diff; data(epoch,track_pair).mean_rate_diff'];
            remapping(epoch,track_pair).mean_max_FR_replay_diff= [remapping(epoch,track_pair).mean_max_FR_replay_diff;  data(epoch,track_pair).mean_max_FR_replay_diff];
            
            % Remapping - track data structure
            remapping_raw(epoch,track_pair).raw_peak_BAYESIAN_plfield_1 = [remapping_raw(epoch,track_pair).raw_peak_BAYESIAN_plfield_1 data(epoch,track_pair).raw_peak_BAYESIAN_plfield_1];
            remapping_raw(epoch,track_pair).raw_peak_BAYESIAN_plfield_2 = [remapping_raw(epoch,track_pair).raw_peak_BAYESIAN_plfield_2 data(epoch,track_pair).raw_peak_BAYESIAN_plfield_2];
            remapping_raw(epoch,track_pair).track1_mean_replay_spikes = [remapping_raw(epoch,track_pair).track1_mean_replay_spikes; data(epoch,track_pair).track1_mean_replay_spikes];
            remapping_raw(epoch,track_pair).track2_mean_replay_spikes = [remapping_raw(epoch,track_pair).track2_mean_replay_spikes; data(epoch,track_pair).track2_mean_replay_spikes];
            remapping_raw(epoch,track_pair).track1_median_replay_rate = [remapping_raw(epoch,track_pair).track1_median_replay_rate; data(epoch,track_pair).track1_median_replay_rate];
            remapping_raw(epoch,track_pair).track2_median_replay_rate=  [remapping_raw(epoch,track_pair).track2_median_replay_rate; data(epoch,track_pair).track2_median_replay_rate];
            remapping_raw(epoch,track_pair).track1_mean_replay_spikes_nonZero = [remapping_raw(epoch,track_pair).track1_mean_replay_spikes_nonZero; data(epoch,track_pair).track1_mean_replay_spikes_nonZero];
            remapping_raw(epoch,track_pair).track2_mean_replay_spikes_nonZero = [remapping_raw(epoch,track_pair).track2_mean_replay_spikes_nonZero; data(epoch,track_pair).track2_mean_replay_spikes_nonZero];
            remapping_raw(epoch,track_pair).track1_median_replay_spikes_nonZero = [remapping_raw(epoch,track_pair).track1_median_replay_spikes_nonZero; data(epoch,track_pair).track1_median_replay_spikes_nonZero];
            remapping_raw(epoch,track_pair).track2_median_replay_spikes_nonZero = [remapping_raw(epoch,track_pair).track2_median_replay_spikes_nonZero; data(epoch,track_pair).track2_median_replay_spikes_nonZero];
            remapping_raw(epoch,track_pair).proportion_events_cell_active_track1 = [remapping_raw(epoch,track_pair).proportion_events_cell_active_track1; data(epoch,track_pair).proportion_events_cell_active_track1];
            remapping_raw(epoch,track_pair).proportion_events_cell_active_track2 = [remapping_raw(epoch,track_pair).proportion_events_cell_active_track2; data(epoch,track_pair).proportion_events_cell_active_track2];
            remapping_raw(epoch,track_pair).mean_rate_T1 = [remapping_raw(epoch,track_pair).mean_rate_T1; data(epoch,track_pair).mean_rate_T1'];
            remapping_raw(epoch,track_pair).mean_rate_T2 = [remapping_raw(epoch,track_pair).mean_rate_T2; data(epoch,track_pair).mean_rate_T2'];
            remapping_raw(epoch,track_pair).proportion_events_T1 = [remapping_raw(epoch,track_pair).proportion_events_T1 data(epoch,track_pair).proportion_events_T1];
            remapping_raw(epoch,track_pair).proportion_events_T2 = [remapping_raw(epoch,track_pair).proportion_events_T2 data(epoch,track_pair).proportion_events_T2];
            remapping_raw(epoch,track_pair).replay_spikes1 = [remapping_raw(epoch,track_pair).replay_spikes1 {data(epoch,track_pair).spikes1}];
            remapping_raw(epoch,track_pair).replay_spikes2 = [remapping_raw(epoch,track_pair).replay_spikes2 {data(epoch,track_pair).spikes2}];
            remapping_raw(epoch,track_pair).replay_rate1 = [remapping_raw(epoch,track_pair).replay_rate1 {data(epoch,track_pair).rate1}];
            remapping_raw(epoch,track_pair).replay_rate2 = [remapping_raw(epoch,track_pair).replay_rate2 {data(epoch,track_pair).rate2}];
            remapping_raw(epoch,track_pair).max_instantaneous_FR_1 = [remapping_raw(epoch,track_pair).max_instantaneous_FR_1 {data(epoch,track_pair).max_instantaneous_FR_1}];
            remapping_raw(epoch,track_pair).max_instantaneous_FR_2 = [remapping_raw(epoch,track_pair).max_instantaneous_FR_2 {data(epoch,track_pair).max_instantaneous_FR_2}];
            remapping_raw(epoch,track_pair).track1_mean_replay_inst_FR_nonZero = [remapping_raw(epoch,track_pair).track1_mean_replay_inst_FR_nonZero; data(epoch,track_pair).track1_mean_replay_inst_FR_nonZero];
            remapping_raw(epoch,track_pair).track2_mean_replay_inst_FR_nonZero = [remapping_raw(epoch,track_pair).track2_mean_replay_inst_FR_nonZero; data(epoch,track_pair).track2_mean_replay_inst_FR_nonZero];

%           remapping(epoch,track_pair).peak_FR_T1=  [remapping(epoch,track_pair).peak_FR_T1; data(epoch,track_pair).peak_FR_T1'];
%           remapping(epoch,track_pair).peak_FR_T2=  [remapping(epoch,track_pair).peak_FR_T2; data(epoch,track_pair).peak_FR_T2'];
             
        end
    end
   cd(master_folder);
end

if ~isempty(folders) && save_option == 1
    switch method
        case 'wcorr'
            save rate_remapping_analysis_TRACK_PAIRS_wcorr remapping remapping_raw
        case 'spearman'
            save rate_remapping_analysis_TRACK_PAIRS_spearman remapping remapping_raw
        case 'control_fixed_spike'
            save rate_remapping_analysis_TRACK_PAIRS_wcorr_FIXED remapping remapping_raw
        case 'rate_detection_control'
            save rate_remapping_analysis_TRACK_PAIRS_wcorr_RATE remapping remapping_raw
        case 'rate_intrinsic_bias_control'
            save rate_remapping_analysis_TRACK_PAIRS_wcorr_INTRINSIC_RATE remapping remapping_raw
        case 'replay_rate_shuffle_control'
            save rate_remapping_analysis_TRACK_PAIRS_wcorr_REPLAY_RATE remapping remapping_raw
        case 'replay_rate_shuffle_detection'
            save rate_remapping_analysis_TRACK_PAIRS_wcorr_REPLAY_RATE_DETECTION remapping remapping_raw
        otherwise
            save rate_remapping_analysis_TRACK_PAIRS remapping remapping_raw
    end
end

end

