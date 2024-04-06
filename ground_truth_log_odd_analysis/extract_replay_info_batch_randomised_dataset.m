function extract_replay_info_batch_randomised_dataset(folders)

% % Extract information cross experiment shuffles
% shuffles={'PRE spike_train_circular_shift','PRE place_field_circular_shift','POST place bin circular shift'...
%     ,'POST time bin permutation','PRE cell_id_shuffle'};
% 
% for n = 1:length(shuffles)
%     [log_odd] = extract_ground_truth_info(folders,'cross_experiment_shuffled','wcorr',shuffles{n},3,2.1)
% 
%     if n == 1
%         save('ground_truth_original\log_odd_wcorr_spike_train_circular_shift_cross_experiment_shuffled.mat','log_odd','-mat')
%     elseif n == 2
%         save('ground_truth_original\log_odd_wcorr_place_field_circular_shift_cross_experiment_shuffled.mat','log_odd','-mat')
%     elseif n == 3
%         save('ground_truth_original\log_odd_wcorr_place_bin_circular_shift_cross_experiment_shuffled.mat','log_odd','-mat')
%     elseif n == 4
%         save('ground_truth_original\log_odd_wcorr_time_bin_permutation_cross_experiment_shuffled.mat','log_odd','-mat')
%     elseif n == 5
%         save('ground_truth_original\log_odd_wcorr_cell_id_shuffle_cross_experiment_shuffled.mat','log_odd','-mat')
%     end
% 
% end
% 
% 
% % place field shuffle log odds
% shuffles={'PRE spike_train_circular_shift','PRE place_field_circular_shift','POST place bin circular shift'...
%     ,'POST time bin permutation','PRE cell_id_shuffle'};
% 
% for n = 1:length(shuffles)
%     [log_odd] = extract_ground_truth_info(folders,'place_field_shifted','wcorr',shuffles{n},3,2.1)
%     if n == 1
%         save('ground_truth_original\log_odd_wcorr_spike_train_circular_shift_place_field_shifted.mat','log_odd','-mat')
%     elseif n == 2
%         save('ground_truth_original\log_odd_wcorr_place_field_circular_shift_place_field_shifted.mat','log_odd','-mat')
%     elseif n == 3
%         save('ground_truth_original\log_odd_wcorr_place_bin_circular_shift_place_field_shifted.mat','log_odd','-mat')
%     elseif n == 4
%         save('ground_truth_original\log_odd_wcorr_time_bin_permutation_place_field_shifted.mat','log_odd','-mat')
%     elseif n == 5
%         save('ground_truth_original\log_odd_wcorr_cell_id_shuffle_place_field_shifted.mat','log_odd','-mat')
%     end
% 
%     clear log_odd
% end
% 
% % spike train shuffle log odds
% shuffles={'PRE spike_train_circular_shift','PRE place_field_circular_shift','POST place bin circular shift'...
%     ,'POST time bin permutation','PRE cell_id_shuffle'};
% 
% for n = 1:length(shuffles)
%     [log_odd] = extract_ground_truth_info(folders,'spike_train_shifted','wcorr',shuffles{n},3,2.1)
% 
%     if n == 1
%         save('ground_truth_original\log_odd_wcorr_spike_train_circular_shift_spike_train_shifted.mat','log_odd','-mat')
%     elseif n == 2
%         save('ground_truth_original\log_odd_wcorr_place_field_circular_shift_spike_train_shifted.mat','log_odd','-mat')
%     elseif n == 3
%         save('ground_truth_original\log_odd_wcorr_place_bin_circular_shift_spike_train_shifted.mat','log_odd','-mat')
%     elseif n == 4
%         save('ground_truth_original\log_odd_wcorr_time_bin_permutation_spike_train_shifted.mat','log_odd','-mat')
%     elseif n == 5
%         save('ground_truth_original\log_odd_wcorr_cell_id_shuffle_spike_train_shifted.mat','log_odd','-mat')
%     end
% end


% Extract place bin circular shift linear and path 
shuffles={'PRE spike_train_circular_shift','PRE place_field_circular_shift','POST place bin circular shift'...
    ,'POST time bin permutation','PRE cell_id_shuffle'};

% cross experiment shuffled randomised dataset
[log_odd] = extract_ground_truth_info(folders,'cross_experiment_shuffled','linear',shuffles{3},3,2.1)
cd ground_truth_original
save log_odd_linear_place_bin_circular_shift_cross_experiment_shuffled log_odd
cd ..

[log_odd] = extract_ground_truth_info(folders,'cross_experiment_shuffled','path',shuffles{3},3,2.1)
cd ground_truth_original
save log_odd_linear_place_bin_circular_shift_cross_experiment_shuffled log_odd
cd ..

% spike train shifted randomised dataset
[log_odd] = extract_ground_truth_info(folders,'spike_train_shifted','linear',shuffles{3},3,2.1)
cd ground_truth_original
save log_odd_linear_place_bin_circular_shift_spike_train_shifted log_odd
cd ..

[log_odd] = extract_ground_truth_info(folders,'spike_train_shifted','path',shuffles{3},3,2.1)
cd ground_truth_original
save log_odd_linear_place_bin_circular_shift_spike_train_shifted log_odd
cd ..

% place field shifted randomised dataset
[log_odd] = extract_ground_truth_info(folders,'place_field_shifted','linear',shuffles{3},3,2.1)
cd ground_truth_original
save log_odd_linear_place_bin_circular_shift_place_field_shifted log_odd
cd ..

[log_odd] = extract_ground_truth_info(folders,'place_field_shifted','path',shuffles{3},3,2.1)
cd ground_truth_original
save log_odd_linear_place_bin_circular_shift_place_field_shifted log_odd
cd ..



% Two shuffles wcorr and linear and path

% cross experiment shuffled randomised dataset
shuffles = 'PRE place and POST time';
method = {'wcorr','linear','path'};
for nmethod = 1:length(method)
    [log_odd] = extract_ground_truth_info(folders,'cross_experiment_shuffled',method{nmethod},shuffles,3,2.1)
    cd ground_truth_original
    if strcmp(method{nmethod},'wcorr')
        save log_odd_wcorr_PRE_place_POST_time_cross_experiment_shuffled log_odd
    elseif strcmp(method{nmethod},'linear')
        save log_odd_linear_PRE_place_POST_time_cross_experiment_shuffled log_odd
    elseif strcmp(method{nmethod},'path')
        save log_odd_path_PRE_place_POST_time_cross_experiment_shuffled log_odd
    end
    cd ..
end

% spike train shifted randomised dataset
shuffles = 'PRE place and POST time';
method = {'wcorr','linear','path'};
for nmethod = 1:length(method)
    [log_odd] = extract_ground_truth_info(folders,'spike_train_shifted',method{nmethod},shuffles,3,2.1)
    cd ground_truth_original
    if strcmp(method{nmethod},'wcorr')
        save log_odd_wcorr_PRE_place_POST_time_spike_train_shifted log_odd
    elseif strcmp(method{nmethod},'linear')
        save log_odd_linear_PRE_place_POST_time_spike_train_shifted log_odd
    elseif strcmp(method{nmethod},'path')
        save log_odd_path_PRE_place_POST_time_spike_train_shifted log_odd
    end
    cd ..
end

% place field shifted randomised dataset
shuffles = 'PRE place and POST time';
method = {'wcorr','linear','path'};
for nmethod = 1:length(method)
    [log_odd] = extract_ground_truth_info(folders,'place_field_shifted',method{nmethod},shuffles,3,2.1)
    cd ground_truth_original
    if strcmp(method{nmethod},'wcorr')
        save log_odd_wcorr_PRE_place_POST_time_place_field_shifted log_odd
    elseif strcmp(method{nmethod},'linear')
        save log_odd_linear_PRE_place_POST_time_place_field_shifted log_odd
    elseif strcmp(method{nmethod},'path')
        save log_odd_path_PRE_place_POST_time_place_field_shifted log_odd
    end
    cd ..
end


% 3 shuffles wcorr, spearman and spearman all spikes
method = {'wcorr','spearman','spearman_all_spikes','linear','path'};
% method = {'wcorr'};
% method = {'linear'};

% cross experiment shuffled randomised dataset
for k = 1:3
    % option = 'global remapped';
    ripple_zscore_threshold = 3;
    no_of_shuffles = 3;
    option = 'cross_experiment_shuffled';
    [log_odd] = extract_ground_truth_info(folders,option,method{k},no_of_shuffles,ripple_zscore_threshold,2.1)

    cd ground_truth_original
    if strcmp(method{k},'wcorr')
        save log_odd_wcorr_cross_experiment_shuffled log_odd
    elseif strcmp(method{k},'spearman')
        save log_odd_spearman_cross_experiment_shuffled log_odd
    elseif strcmp(method{k},'spearman_all_spikes')
        save log_odd_spearman_all_spikes_cross_experiment_shuffled log_odd
%     elseif strcmp(method{k},'linear')
%         save log_odd_linear_cross_experiment_shuffled log_odd
%     elseif strcmp(method{k},'path')
%         save log_odd_path_cross_experiment_shuffled log_odd
    end
    cd ..
end


% spike train shifted randomised dataset
for k = 1:3
    % option = 'global remapped';
    ripple_zscore_threshold = 3;
    no_of_shuffles = 3;
    option = 'spike_train_shifted';
    [log_odd] = extract_ground_truth_info(folders,option,method{k},no_of_shuffles,ripple_zscore_threshold,2.1)

    cd ground_truth_original
    if strcmp(method{k},'wcorr')
        save log_odd_wcorr_spike_train_shifted log_odd
    elseif strcmp(method{k},'spearman')
        save log_odd_spearman_spike_train_shifted log_odd
    elseif strcmp(method{k},'spearman_all_spikes')
        save log_odd_spearman_all_spikes_spike_train_shifted log_odd
%     elseif strcmp(method{k},'linear')
%         save log_odd_linear_cross_experiment_shuffled log_odd
%     elseif strcmp(method{k},'path')
%         save log_odd_path_cross_experiment_shuffled log_odd
    end
    cd ..
end


% place field shifted randomised dataset
for k = 1:3
    % option = 'global remapped';
    ripple_zscore_threshold = 3;
    no_of_shuffles = 3;
    option = 'place_field_shifted';
    [log_odd] = extract_ground_truth_info(folders,option,method{k},no_of_shuffles,ripple_zscore_threshold,2.1)

    cd ground_truth_original
    if strcmp(method{k},'wcorr')
        save log_odd_wcorr_place_field_shifted log_odd
    elseif strcmp(method{k},'spearman')
        save log_odd_spearman_place_field_shifted log_odd
    elseif strcmp(method{k},'spearman_all_spikes')
        save log_odd_spearman_all_spikes_place_field_shifted log_odd
%     elseif strcmp(method{k},'linear')
%         save log_odd_linear_cross_experiment_shuffled log_odd
%     elseif strcmp(method{k},'path')
%         save log_odd_path_cross_experiment_shuffled log_odd
    end
    cd ..
end

end