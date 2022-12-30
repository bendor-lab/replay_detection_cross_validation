function extract_replay_info_batch(folders)

shuffles={'PRE spike_train_circular_shift','PRE place_field_circular_shift','POST place bin circular shift','POST time bin circular shift','POST time bin permutation'};

% Spike shuffle
[log_odd] = extract_ground_truth_info(folders,'original','wcorr',shuffles{1},3,2.1)
cd ground_truth_original
save log_odd_wcorr_spike_train_circular_shift log_odd
cd ..

[log_odd] = extract_ground_truth_info(folders,'global remapped','wcorr',shuffles{1},3,2.1)
cd ground_truth_original
save log_odd_wcorr_spike_train_circular_shift_global_remapped log_odd
cd ..

% Place field shuffle
[log_odd] = extract_ground_truth_info(folders,'original','wcorr',shuffles{2},3,2.1)
cd ground_truth_original
save log_odd_wcorr_place_field_circular_shift log_odd
cd ..

[log_odd] = extract_ground_truth_info(folders,'global remapped','wcorr',shuffles{2},3,2.1)
cd ground_truth_original
save log_odd_wcorr_place_field_circular_shift_global_remapped log_odd
cd ..

% Place bin shuffle (wcorr)
[log_odd] = extract_ground_truth_info(folders,'original','wcorr',shuffles{3},3,2.1)
cd ground_truth_original
save log_odd_wcorr_place_bin_circular_shift log_odd
cd ..

[log_odd] = extract_ground_truth_info(folders,'global remapped','wcorr',shuffles{3},3,2.1)
cd ground_truth_original
save log_odd_wcorr_place_bin_circular_shift_global_remapped log_odd
cd ..

% Time bin shuffle
[log_odd] = extract_ground_truth_info(folders,'original','wcorr',shuffles{5},3,2.1)
cd ground_truth_original
save log_odd_wcorr_time_bin_permutation log_odd
cd ..

[log_odd] = extract_ground_truth_info(folders,'global remapped','wcorr',shuffles{5},3,2.1)
cd ground_truth_original
save log_odd_wcorr_time_bin_permutation_global_remapped log_odd
cd ..


% Place bin shuffle (linear)
[log_odd] = extract_ground_truth_info(folders,'original','linear',shuffles{3},3,2.1)
cd ground_truth_original
save log_odd_linear_place_bin_circular_shift log_odd
cd ..

[log_odd] = extract_ground_truth_info(folders,'global remapped','linear',shuffles{3},3,2.1)
cd ground_truth_original
save log_odd_linear_place_bin_circular_shift_global_remapped log_odd
cd ..

% Place bin shuffle (path)
[log_odd] = extract_ground_truth_info(folders,'original','path',shuffles{3},3,2.1)
cd ground_truth_original
save log_odd_path_place_bin_circular_shift log_odd
cd ..

[log_odd] = extract_ground_truth_info(folders,'global remapped','linear',shuffles{3},3,2.1)
cd ground_truth_original
save log_odd_path_place_bin_circular_shift_global_remapped log_odd
cd ..


% Shuffle 1 + 2
shuffles = 'two PRE';
[log_odd] = extract_ground_truth_info(folders,'original','wcorr',shuffles,3,2.1)
cd ground_truth_original
save log_odd_wcorr_two_PRE_shuffles log_odd
cd ..

[log_odd] = extract_ground_truth_info(folders,'global remapped','wcorr',shuffles,3,2.1)
cd ground_truth_original
save log_odd_wcorr_two_PRE_shuffles_global_remapped log_odd
cd ..

% Shuffle 3 + 5 
shuffles = 'two POST'
[log_odd] = extract_ground_truth_info(folders,'original','wcorr',shuffles,3,2.1)
cd ground_truth_original
save log_odd_wcorr_two_POST_shuffles log_odd
cd ..

[log_odd] = extract_ground_truth_info(folders,'global remapped','wcorr',shuffles,3,2.1)
cd ground_truth_original
save log_odd_wcorr_two_POST_shuffles_global_remapped log_odd
cd ..

% Shuffle 2 + 5
shuffles = 'PRE place and POST time';
method = {'wcorr','linear','path'};

for nmethod = 1:length(method)
    [log_odd] = extract_ground_truth_info(folders,'original',method{nmethod},shuffles,3,2.1)
    cd ground_truth_original
    if strcmp(method{nmethod},'wcorr')
        save log_odd_wcorr_PRE_place_POST_time log_odd
    elseif strcmp(method{nmethod},'linear')
        save log_odd_linear_PRE_place_POST_time log_odd
    elseif strcmp(method{nmethod},'path')
        save log_odd_path_PRE_place_POST_time log_odd
    end
    cd ..
end


method = {'wcorr','linear','path'};
for nmethod = 1:length(method)
    [log_odd] = extract_ground_truth_info(folders,'global remapped',method{nmethod},shuffles,3,2.1)
    cd ground_truth_original
    if strcmp(method{nmethod},'wcorr')
        save log_odd_wcorr_PRE_place_POST_time_global_remapped log_odd
    elseif strcmp(method{nmethod},'linear')
        save log_odd_linear_PRE_place_POST_time_global_remapped log_odd
    elseif strcmp(method{nmethod},'path')
        save log_odd_path_PRE_place_POST_time_global_remapped log_odd
    end
    cd ..
end


% % Shuffle 2 + 3
shuffles = 'PRE place and POST place';
[log_odd] = extract_ground_truth_info(folders,'original','wcorr',shuffles,3,2.1)
cd ground_truth_original
save log_odd_wcorr_PRE_place_POST_place log_odd
cd ..

[log_odd] = extract_ground_truth_info(folders,'global remapped','wcorr',shuffles,3,2.1)
cd ground_truth_original
save log_odd_wcorr_PRE_place_POST_place_global_remapped log_odd
cd ..


% Shuffle 1 + 3
shuffles = 'PRE spike and POST place';
[log_odd] = extract_ground_truth_info(folders,'original','wcorr',shuffles,3,2.1)
cd ground_truth_original
save log_odd_wcorr_PRE_spike_POST_place log_odd
cd ..


[log_odd] = extract_ground_truth_info(folders,'global remapped','wcorr',shuffles,3,2.1)
cd ground_truth_original
save log_odd_wcorr_PRE_spike_POST_place_global_remapped log_odd
cd ..


% Ripple threshold 0 
ripple_zscore_threshold = 0;
no_of_shuffles = 3;
[log_odd] = extract_ground_truth_info(folders,'original','wcorr',no_of_shuffles,ripple_zscore_threshold,2.01)
cd ground_truth_original
save log_odd_wcorr_ripple_0 log_odd
cd ..

ripple_zscore_threshold = 0;
no_of_shuffles = 3;
[log_odd] = extract_ground_truth_info(folders,'global remapped','wcorr',no_of_shuffles,ripple_zscore_threshold,2.01)
cd ground_truth_original
save log_odd_wcorr_ripple_0_global_remapped log_odd
cd ..

% Ripple threshold 0 - 2 shuffles
shuffles = 'PRE place and POST time';
ripple_zscore_threshold = 0;
[log_odd] = extract_ground_truth_info(folders,'original','wcorr',shuffles,ripple_zscore_threshold,2.01)
cd ground_truth_original
save log_odd_wcorr_PRE_place_POST_time_ripple_0 log_odd
cd ..

ripple_zscore_threshold = 0;
[log_odd] = extract_ground_truth_info(folders,'global remapped','wcorr',shuffles,ripple_zscore_threshold,2.01)
cd ground_truth_original
save log_odd_wcorr_PRE_place_POST_time_ripple_0_global_remapped log_odd
cd ..

% Ripple threshold 0 - 1 shuffle
shuffles = 'POST place bin circular shift';
ripple_zscore_threshold = 0;
[log_odd] = extract_ground_truth_info(folders,'original','wcorr',shuffles,ripple_zscore_threshold,2.01)
cd ground_truth_original
save log_odd_wcorr_POST_place_ripple_0 log_odd
cd ..

ripple_zscore_threshold = 0;
[log_odd] = extract_ground_truth_info(folders,'global remapped','wcorr',shuffles,ripple_zscore_threshold,2.01)
cd ground_truth_original
save log_odd_wcorr_POST_place_ripple_0_global_remapped log_odd
cd ..

method = {'wcorr','spearman','spearman_all_spikes','linear','path'};
% method = {'wcorr'};
% method = {'linear'};
% Original ripple threshold 3
for k = 1:length(method)
    % option = 'global remapped';
    ripple_zscore_threshold = 3;
    option = 'original';
    
    no_of_shuffles = 3;
    [log_odd] = extract_ground_truth_info(folders,option,method{k},no_of_shuffles,ripple_zscore_threshold,2.1)

    cd ground_truth_original
    if strcmp(method{k},'wcorr')
        save log_odd_wcorr log_odd
    elseif strcmp(method{k},'spearman')
        save log_odd_spearman log_odd
    elseif strcmp(method{k},'spearman_all_spikes')
        save log_odd_spearman_all_spikes log_odd
    elseif strcmp(method{k},'linear')
        save log_odd_linear log_odd
    elseif strcmp(method{k},'path')
        save log_odd_path log_odd
    end
    cd ..


    no_of_shuffles = 3;
    option = 'global remapped';
    [log_odd] = extract_ground_truth_info(folders,option,method{k},no_of_shuffles,ripple_zscore_threshold,2.1)

    cd ground_truth_original
    if strcmp(method{k},'wcorr')
        save log_odd_wcorr_global_remapped log_odd
    elseif strcmp(method{k},'spearman')
        save log_odd_spearman_global_remapped log_odd
    elseif strcmp(method{k},'spearman_all_spikes')
        save log_odd_spearman_all_spikes_global_remapped log_odd
    elseif strcmp(method{k},'linear')
        save log_odd_linear_global_remapped log_odd
    elseif strcmp(method{k},'path')
        save log_odd_path_global_remapped log_odd
    end
    cd ..
end



end