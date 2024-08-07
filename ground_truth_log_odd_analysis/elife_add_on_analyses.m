%% Code for additional analyses

%% creating and analysing randomised dataset 
clear all
workingDir = 'P:\ground_truth_replay_analysis\Dropbo_data8';
cd(workingDir)
folders = {'2019-06-20_09-55-42' 'N-BLU_Day7_Ctrl-16x30_no_reexp' 'N-BLU_Day8_Ctrl-15-NoRest-15' 'Q-BLU_Day2_RateRemap' 'Q-BLU_Day3_RateRemap'...
    'RAT1_2018-10-05_09-42-15' 'RAT4_2019-06-15_10-46-59' 'RAT4_2019-06-19_10-00-07' 'RAT5_2019-06-24_10-43-24' 'RAT5_2019-06-26_10-27-15'};

BAYSESIAN_NORMALIZED_ACROSS_TRACKS = 1; % Normalised across tracks
timebin = 0.02;
% ground_truth_log_odd_analysis_place_field_randomised(folders,0.02,BAYSESIAN_NORMALIZED_ACROSS_TRACKS)
ground_truth_log_odd_analysis_spike_train_randomised(folders,0.02,BAYSESIAN_NORMALIZED_ACROSS_TRACKS)
% ground_truth_log_odd_analysis_spike_train_randomised_addon(folders,0.02,BAYSESIAN_NORMALIZED_ACROSS_TRACKS)
% place_fields_BAYESIAN_combined = generate_cross_experiment_cell_id_shuffles(folders);
% ground_truth_log_odd_analysis_cross_experiment_randomised(folders,0.02,BAYSESIAN_NORMALIZED_ACROSS_TRACKS)

folders = {'2019-06-20_09-55-42' 'N-BLU_Day7_Ctrl-16x30_no_reexp' 'N-BLU_Day8_Ctrl-15-NoRest-15' 'Q-BLU_Day2_RateRemap' 'Q-BLU_Day3_RateRemap'...
    'RAT1_2018-10-05_09-42-15' 'RAT4_2019-06-15_10-46-59' 'RAT4_2019-06-19_10-00-07' 'RAT5_2019-06-24_10-43-24' 'RAT5_2019-06-26_10-27-15'};

%% Extract information and save structures
clear all
workingDir = 'P:\ground_truth_replay_analysis\Dropbo_data8';
cd(workingDir)
folders = {'2019-06-20_09-55-42' 'N-BLU_Day7_Ctrl-16x30_no_reexp' 'N-BLU_Day8_Ctrl-15-NoRest-15' 'Q-BLU_Day2_RateRemap' 'Q-BLU_Day3_RateRemap'...
    'RAT1_2018-10-05_09-42-15' 'RAT4_2019-06-15_10-46-59' 'RAT4_2019-06-19_10-00-07' 'RAT5_2019-06-24_10-43-24' 'RAT5_2019-06-26_10-27-15'};

calculate_jump_distance(folders,'place_field_shifted')
calculate_jump_distance(folders,'spike_train_shifted')
calculate_jump_distance(folders,'cross_experiment_shuffled')


% Re running spearman corr for randomised dataset
workingDir = 'D:\Dropbo_data8';
spearman_shuffled_data_re_analysis(folders,'cross_experiment_shuffled');
spearman_shuffled_data_re_analysis(folders,'spike_train_shifted');
cd(workingDir)
spearman_shuffled_data_re_analysis(folders,'place_field_shifted');

% shuffles={'PRE spike_train_circular_shift','PRE place_field_circular_shift','POST place bin circular shift'...
%     ,'POST time bin permutation','PRE cell_id_shuffle'};


clear all
workingDir = 'D:\Dropbo_data8';
cd(workingDir)
folders = {'2019-06-20_09-55-42' 'N-BLU_Day7_Ctrl-16x30_no_reexp' 'N-BLU_Day8_Ctrl-15-NoRest-15' 'Q-BLU_Day2_RateRemap' 'Q-BLU_Day3_RateRemap'...
    'RAT1_2018-10-05_09-42-15' 'RAT4_2019-06-15_10-46-59' 'RAT4_2019-06-19_10-00-07' 'RAT5_2019-06-24_10-43-24' 'RAT5_2019-06-26_10-27-15'};

% extract_replay_info_batch_cross_experiment_shuffled(folders)
extract_replay_info_batch_randomised_dataset(folders)

%% Validating the use of cell-id randomized dataset
clear all
% workingDir = 'D:\ground_truth_replay_analysis\Dropbo_data8_normalised_within';
workingDir = 'P:\ground_truth_replay_analysis\Dropbo_data8';
cd(workingDir)
folders = {'2019-06-20_09-55-42' 'N-BLU_Day7_Ctrl-16x30_no_reexp' 'N-BLU_Day8_Ctrl-15-NoRest-15' 'Q-BLU_Day2_RateRemap' 'Q-BLU_Day3_RateRemap'...
    'RAT1_2018-10-05_09-42-15' 'RAT4_2019-06-15_10-46-59' 'RAT4_2019-06-19_10-00-07' 'RAT5_2019-06-24_10-43-24' 'RAT5_2019-06-26_10-27-15'};
option = 'common';

method = {'spike shuffle','place field shuffle','place shuffle','time shuffle','cell id shuffle'};
plot_ground_truth_cell_id_shuffle_validation(folders,option,method)

%% dual-track vs false positive rates
clear all
workingDir = 'P:\ground_truth_replay_analysis\Dropbo_data8';
option = 'common';
cd(workingDir)
folders = {'2019-06-20_09-55-42' 'N-BLU_Day7_Ctrl-16x30_no_reexp' 'N-BLU_Day8_Ctrl-15-NoRest-15' 'Q-BLU_Day2_RateRemap' 'Q-BLU_Day3_RateRemap'...
    'RAT1_2018-10-05_09-42-15' 'RAT4_2019-06-15_10-46-59' 'RAT4_2019-06-19_10-00-07' 'RAT5_2019-06-24_10-43-24' 'RAT5_2019-06-26_10-27-15'};
method = {'wcorr 1 shuffle','wcorr 1 shuffle + jump distance','wcorr 2 shuffles','wcorr 3 shuffles'...
    'spearman median spike','spearman all spikes',...
    'linear 1 shuffle','linear 2 shuffles'};
plot_ground_truth_compare_multitrack(folders,option,method)

%% Compare Cell ID randomisation vs place field randomisation vs spike train
% randomisation
workingDir = 'D:\Dropbo_data8';
cd(workingDir)
folders = {'2019-06-20_09-55-42' 'N-BLU_Day7_Ctrl-16x30_no_reexp' 'N-BLU_Day8_Ctrl-15-NoRest-15' 'Q-BLU_Day2_RateRemap' 'Q-BLU_Day3_RateRemap'...
    'RAT1_2018-10-05_09-42-15' 'RAT4_2019-06-15_10-46-59' 'RAT4_2019-06-19_10-00-07' 'RAT5_2019-06-24_10-43-24' 'RAT5_2019-06-26_10-27-15'};
option = 'common';
% method = {'wcorr','spearman','linear'};
method = {'wcorr 1 shuffle','wcorr 1 shuffle + jump distance','wcorr 2 shuffles','wcorr 3 shuffles'...
    'spearman median spike','spearman all spikes',...
    'linear 1 shuffle','linear 2 shuffles'};
plot_ground_truth_compare_shuffles_with_different_randomisation(folders,option,method)

%% Plotting log odds difference within each session
clear all
% workingDir = 'D:\ground_truth_replay_analysis\Dropbo_data8_normalised_within';
workingDir = 'P:\ground_truth_replay_analysis\Dropbo_data8';
cd(workingDir)

% folders = { '2019-06-20_09-55-42' 'N-BLU_Day7_Ctrl-16x30_no_reexp' 'N-BLU_Day8_Ctrl-15-NoRest-15' 'Q-BLU_Day2_RateRemap' 'Q-BLU_Day3_RateRemap'};
folders = {'2019-06-20_09-55-42' 'N-BLU_Day7_Ctrl-16x30_no_reexp' 'N-BLU_Day8_Ctrl-15-NoRest-15' 'Q-BLU_Day2_RateRemap' 'Q-BLU_Day3_RateRemap'...
    'RAT1_2018-10-05_09-42-15' 'RAT4_2019-06-15_10-46-59' 'RAT4_2019-06-19_10-00-07' 'RAT5_2019-06-24_10-43-24' 'RAT5_2019-06-26_10-27-15'};
option = 'common';
method = {'wcorr 2 shuffles','spearman median spike'};
plot_ground_truth_across_sessions(folders,option,method)

%% Ripple Threshold comparision

clear all
% workingDir = 'D:\ground_truth_replay_analysis\Dropbo_data8_normalised_within';
workingDir = 'P:\ground_truth_replay_analysis\Dropbo_data8';
% workingDir = 'P:\ground_truth_replay_analysis\Dropbo_data8';
cd(workingDir)

% cd('P:\ground_truth_replay_analysis\Dropbo_data8 - copy');
% folders = { '2019-06-20_09-55-42' 'N-BLU_Day7_Ctrl-16x30_no_reexp' 'N-BLU_Day8_Ctrl-15-NoRest-15' 'Q-BLU_Day2_RateRemap' 'Q-BLU_Day3_RateRemap'};
folders = {'2019-06-20_09-55-42' 'N-BLU_Day7_Ctrl-16x30_no_reexp' 'N-BLU_Day8_Ctrl-15-NoRest-15' 'Q-BLU_Day2_RateRemap' 'Q-BLU_Day3_RateRemap'...
    'RAT1_2018-10-05_09-42-15' 'RAT4_2019-06-15_10-46-59' 'RAT4_2019-06-19_10-00-07' 'RAT5_2019-06-24_10-43-24' 'RAT5_2019-06-26_10-27-15'};
option = 'common';

method = {'ripple 0', 'ripple 3','ripple 5','ripple 10'};
plot_ground_truth_compare_ripple_threshold(folders,option,method)

%% Ripple Threshold comparision with lower speed threshold

clear all
% workingDir = 'D:\ground_truth_replay_analysis\Dropbo_data8_normalised_within';
workingDir = 'P:\ground_truth_replay_analysis\Dropbo_data8';
% workingDir = 'P:\ground_truth_replay_analysis\Dropbo_data8';
cd(workingDir)

% cd('P:\ground_truth_replay_analysis\Dropbo_data8 - copy');
% folders = { '2019-06-20_09-55-42' 'N-BLU_Day7_Ctrl-16x30_no_reexp' 'N-BLU_Day8_Ctrl-15-NoRest-15' 'Q-BLU_Day2_RateRemap' 'Q-BLU_Day3_RateRemap'};
folders = {'2019-06-20_09-55-42' 'N-BLU_Day7_Ctrl-16x30_no_reexp' 'N-BLU_Day8_Ctrl-15-NoRest-15' 'Q-BLU_Day2_RateRemap' 'Q-BLU_Day3_RateRemap'...
    'RAT1_2018-10-05_09-42-15' 'RAT4_2019-06-15_10-46-59' 'RAT4_2019-06-19_10-00-07' 'RAT5_2019-06-24_10-43-24' 'RAT5_2019-06-26_10-27-15'};
option = 'common';
method = {'ripple 0-3', 'ripple 3-5','ripple 5-10','ripple 10 and above'};
plot_ground_truth_compare_ripple_replay_low_speed(folders,option,method)


%% Plot 4 individual examples of replay events
clear all
folders = {'2019-06-20_09-55-42' 'N-BLU_Day7_Ctrl-16x30_no_reexp' 'N-BLU_Day8_Ctrl-15-NoRest-15' 'Q-BLU_Day2_RateRemap' 'Q-BLU_Day3_RateRemap'...
    'RAT1_2018-10-05_09-42-15' 'RAT4_2019-06-15_10-46-59' 'RAT4_2019-06-19_10-00-07' 'RAT5_2019-06-24_10-43-24' 'RAT5_2019-06-26_10-27-15'};
workingDir = 'P:\ground_truth_replay_analysis\Dropbo_data8';
cd(fullfile(workingDir,'ground_truth_original/'))
load('log_odd_wcorr_PRE_place_POST_time.mat')

cd(fullfile(workingDir,folders{7}))
load decoded_replay_events.mat

close all
good_events = [3265 1630 4056 3737];
good_sequence_pvalue = [];
good_sequence_log_odds = [];

for event = 1:4
    index(event) = find(log_odd.index == good_events(event) & log_odd.experiment==7);
end

good_sequence_pvalue(1,:) = log_odd.pvalue(1,index);
good_sequence_pvalue(2,:) = log_odd.pvalue(2,index);

good_sequence_log_odds(1,:) = log_odd.common_zscored.original(1,index);
good_sequence_log_odds(2,:) = log_odd.common_zscored.original(2,index);


nfigure = 1;
fig = figure(nfigure);
fig.Position = [435 105 1280 885];
c = 1;
for event = 1:length(good_events)
    if c >= 16
        fig = figure(nfigure);
        fig.Position = [435 105 870 885];
        nfigure = nfigure + 1;
        c = 1;
    end

    
    subplot(4,4,c)
    imagesc(decoded_replay_events(1).replay_events(good_events(event)).decoded_position)
    colormap(flipud(bone))
    colorbar
    caxis([0 max(max([decoded_replay_events(1).replay_events(good_events(event)).decoded_position...
        decoded_replay_events(2).replay_events(good_events(event)).decoded_position]))])
    %     box off
    set(gca,'XTick',[],'YTick',[],"TickDir","out",'Color','none')
    title(sprintf('Session 7 Event %i T1 \nP value = %.2f log odds %.2f',good_events(event),good_sequence_pvalue(1,event),good_sequence_log_odds(1,event)));
    c = c + 1;

    subplot(4,4,c)
    imagesc(decoded_replay_events(2).replay_events(good_events(event)).decoded_position)
    colormap(flipud(bone))
    colorbar
    caxis([0 max(max([decoded_replay_events(1).replay_events(good_events(event)).decoded_position...
        decoded_replay_events(2).replay_events(good_events(event)).decoded_position]))])
    %     box off
    set(gca,'XTick',[],'YTick',[],"TickDir","out",'Color','none')
    %     title(sprintf('Session 7 Event %i T2 \nP value = %.2f log odds %.2f',good_events(event),good_sequence_pvalue(2,event),good_sequence_log_odds(1,event)))
    c = c + 1;
end

cd(fullfile(workingDir,'ground_truth_original/resubmission_figure'))
% cd ground_truth_original\resubmission_figure
filename = sprintf('4 example decoded probability.pdf')
saveas(gcf,filename)
filename = sprintf('4 example decoded probability.fig')
saveas(gcf,filename)
cd(workingDir)
clf
