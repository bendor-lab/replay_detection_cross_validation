function extract_replay_info_batch_cross_experiment_shuffled(folders)

shuffles ={'PRE spike_train_circular_shift','PRE place_field_circular_shift','POST place bin circular shift'...
    ,'POST time bin permutation','PRE cell_id_shuffle'};

for n = 1:length(shuffles)
    [log_odd] = extract_ground_truth_info(folders,'cross_experiment_shuffled','wcorr',shuffles{n},3,2.1)

    if n == 1
        save('ground_truth_original\log_odd_wcorr_spike_train_circular_shift_cross_experiment_shuffled.mat','log_odd','-mat')
    elseif n == 2
        save('ground_truth_original\log_odd_wcorr_place_field_circular_shift_cross_experiment_shuffled.mat','log_odd','-mat')
    elseif n == 3
        save('ground_truth_original\log_odd_wcorr_place_bin_circular_shift_cross_experiment_shuffled.mat','log_odd','-mat')
    elseif n == 4
        save('ground_truth_original\log_odd_wcorr_time_bin_permutation_cross_experiment_shuffled.mat','log_odd','-mat')
    elseif n == 5
        save('ground_truth_original\log_odd_wcorr_cell_id_shuffle_cross_experiment_shuffled.mat','log_odd','-mat')
    end

end

[log_odd] = extract_ground_truth_info(folders,'cross_experiment_shuffled','linear',shuffles{3},3,2.1)
cd ground_truth_original
save log_odd_linear_place_bin_circular_shift_cross_experiment_shuffled log_odd
cd ..


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


method = {'wcorr','spearman','spearman_all_spikes','linear','path'};
% method = {'wcorr'};
% method = {'linear'};
% Original ripple threshold 3
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
    elseif strcmp(method{k},'linear')
        save log_odd_linear_cross_experiment_shuffled log_odd
    elseif strcmp(method{k},'path')
        save log_odd_path_cross_experiment_shuffled log_odd
    end
    cd ..
end


% place field shifted randomised dataset

shuffles ={'PRE spike_train_circular_shift','PRE place_field_circular_shift','POST place bin circular shift'...
    ,'POST time bin permutation','PRE cell_id_shuffle'};
[log_odd] = extract_ground_truth_info(folders,'place_field_shifted','linear',shuffles{3},3,2.1)
cd ground_truth_original
save log_odd_linear_place_bin_circular_shift_place_field_shifted log_odd
cd ..

method = {'wcorr','spearman','spearman_all_spikes','linear','path'};
% method = {'wcorr'};
% method = {'linear'};
% Original ripple threshold 3
for k = 2:3
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
    elseif strcmp(method{k},'linear')
        save log_odd_linear_place_field_shifted log_odd
    elseif strcmp(method{k},'path')
        save log_odd_path_place_field_shifted log_odd
    end
    cd ..
end

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

end