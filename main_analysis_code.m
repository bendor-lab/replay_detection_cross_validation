%% Cross-validation of replay detection performance by using combining both sequenceness and seuqenceless decoding
% Sequenceness analysis

clear all

% workingDir = 'P:\ground_truth_replay_analysis\Dropbo_data8_normalised_within';
workingDir = 'P:\ground_truth_replay_analysis\Dropbo_data8';

% cd P:\ground_truth_replay_analysis\Dropbo_data8_normalised_within\;
cd(workingDir)
% folders = { '2019-06-20_09-55-42' 'N-BLU_Day7_Ctrl-16x30_no_reexp' 'N-BLU_Day8_Ctrl-15-NoRest-15' 'Q-BLU_Day2_RateRemap' 'Q-BLU_Day3_RateRemap'}
folders = {'2019-06-20_09-55-42' 'N-BLU_Day7_Ctrl-16x30_no_reexp' 'N-BLU_Day8_Ctrl-15-NoRest-15' 'Q-BLU_Day2_RateRemap' 'Q-BLU_Day3_RateRemap'...
    'RAT1_2018-10-05_09-42-15' 'RAT4_2019-06-15_10-46-59' 'RAT4_2019-06-19_10-00-07' 'RAT5_2019-06-24_10-43-24' 'RAT5_2019-06-26_10-27-15'};

% folders = {'RAT4_2019-06-15_10-46-59' 'RAT4_2019-06-19_10-00-07' 'RAT5_2019-06-24_10-43-24' 'RAT5_2019-06-26_10-27-15'};

BAYSESIAN_NORMALIZED_ACROSS_TRACKS = 0 % Normalised within tracks
% BAYSESIAN_NORMALIZED_ACROSS_TRACKS = 1 % Normalised across tracks

ground_truth_replay_sequence_analysis(folders,BAYSESIAN_NORMALIZED_ACROSS_TRACKS)
% ground_truth_replay_sequence_analysis_addon(folders,BAYSESIAN_NORMALIZED_ACROSS_TRACKS)
% ground_truth_replay_sequence_analysis_addon_global_remapped(folders,BAYSESIAN_NORMALIZED_ACROSS_TRACKS)

% Sequenceless decoding + global reampped (sequenceness + sequenceless)
clear all
% cd P:\ground_truth_replay_analysis\Dropbo_data8;
cd P:\ground_truth_replay_analysis\Dropbo_data8_normalised_within;

folders = {'2019-06-20_09-55-42' 'N-BLU_Day7_Ctrl-16x30_no_reexp' 'N-BLU_Day8_Ctrl-15-NoRest-15' 'Q-BLU_Day2_RateRemap' 'Q-BLU_Day3_RateRemap'...
    'RAT1_2018-10-05_09-42-15' 'RAT4_2019-06-15_10-46-59' 'RAT4_2019-06-19_10-00-07' 'RAT5_2019-06-24_10-43-24' 'RAT5_2019-06-26_10-27-15'};


% Bayesian Decoding of Track 1 and Track 2 Replay -> Log Odd
% BAYSESIAN_NORMALIZED_ACROSS_TRACKS = 0 % Normalised within tracks
BAYSESIAN_NORMALIZED_ACROSS_TRACKS = 1 % Normalised across tracks
ground_truth_log_odd_analysis(folders,0.02,BAYSESIAN_NORMALIZED_ACROSS_TRACKS);

if ~isfolder({'ground_truth_original'})
    mkdir('ground_truth_original')
end



%% Re-analysing spearman for all spikes and cell id shuffled
clear all
workingDir = 'P:\ground_truth_replay_analysis\Dropbo_data8_normalised_within';
% workingDir = 'P:\ground_truth_replay_analysis\Dropbo_data8';

cd(workingDir)

folders = {'2019-06-20_09-55-42' 'N-BLU_Day7_Ctrl-16x30_no_reexp' 'N-BLU_Day8_Ctrl-15-NoRest-15' 'Q-BLU_Day2_RateRemap' 'Q-BLU_Day3_RateRemap'...
    'RAT1_2018-10-05_09-42-15' 'RAT4_2019-06-15_10-46-59' 'RAT4_2019-06-19_10-00-07' 'RAT5_2019-06-24_10-43-24' 'RAT5_2019-06-26_10-27-15'};
spearman_analysis(folders)

folders = {'2019-06-20_09-55-42' 'N-BLU_Day7_Ctrl-16x30_no_reexp' 'N-BLU_Day8_Ctrl-15-NoRest-15' 'Q-BLU_Day2_RateRemap' 'Q-BLU_Day3_RateRemap'...
    'RAT1_2018-10-05_09-42-15' 'RAT4_2019-06-15_10-46-59' 'RAT4_2019-06-19_10-00-07' 'RAT5_2019-06-24_10-43-24' 'RAT5_2019-06-26_10-27-15'};

spearman_shuffled_data_re_analysis(folders,'global_remapped')

method = {'wcorr','spearman','linear','path'};
folders = {'2019-06-20_09-55-42' 'N-BLU_Day7_Ctrl-16x30_no_reexp' 'N-BLU_Day8_Ctrl-15-NoRest-15' 'Q-BLU_Day2_RateRemap' 'Q-BLU_Day3_RateRemap'...
    'RAT1_2018-10-05_09-42-15' 'RAT4_2019-06-15_10-46-59' 'RAT4_2019-06-19_10-00-07' 'RAT5_2019-06-24_10-43-24' 'RAT5_2019-06-26_10-27-15'};

cd(workingDir)
calculate_jump_distance(folders,'common')
cd(workingDir)
calculate_jump_distance(folders,'global_remapped')


clear all
workingDir = 'D:\ground_truth_replay_analysis\Dropbo_data8_normalised_within';
% workingDir = 'P:\ground_truth_replay_analysis\Dropbo_data8';
cd(workingDir)
folders = {'2019-06-20_09-55-42' 'N-BLU_Day7_Ctrl-16x30_no_reexp' 'N-BLU_Day8_Ctrl-15-NoRest-15' 'Q-BLU_Day2_RateRemap' 'Q-BLU_Day3_RateRemap'...
    'RAT1_2018-10-05_09-42-15' 'RAT4_2019-06-15_10-46-59' 'RAT4_2019-06-19_10-00-07' 'RAT5_2019-06-24_10-43-24' 'RAT5_2019-06-26_10-27-15'};

extract_replay_info_batch(folders);

% for n = 1:length(folders)
%     cd(folders{n})
%     folders{n}
%     load global_remapped_common_good_probability_ratio
%     load global_remapped_original_probability_ratio
%     cd ..
% end

%% Calculate jump distance and Extract log odd information
clear all

% workingDir = 'P:\ground_truth_replay_analysis\Dropbo_data8_normalised_within';
workingDir = 'P:\ground_truth_replay_analysis\Dropbo_data8';
cd(workingDir)

method = {'wcorr','spearman','linear','path'};
folders = {'2019-06-20_09-55-42' 'N-BLU_Day7_Ctrl-16x30_no_reexp' 'N-BLU_Day8_Ctrl-15-NoRest-15' 'Q-BLU_Day2_RateRemap' 'Q-BLU_Day3_RateRemap'...
    'RAT1_2018-10-05_09-42-15' 'RAT4_2019-06-15_10-46-59' 'RAT4_2019-06-19_10-00-07' 'RAT5_2019-06-24_10-43-24' 'RAT5_2019-06-26_10-27-15'};

calculate_jump_distance(folders,'common')
calculate_jump_distance(folders,'global remapped')
extract_replay_info_batch(folders);


clear all
workingDir = 'D:\ground_truth_replay_analysis\Dropbo_data8_normalised_within';
% workingDir = 'P:\ground_truth_replay_analysis\Dropbo_data8';
cd(workingDir)

method = {'wcorr','spearman','linear','path'};
folders = {'2019-06-20_09-55-42' 'N-BLU_Day7_Ctrl-16x30_no_reexp' 'N-BLU_Day8_Ctrl-15-NoRest-15' 'Q-BLU_Day2_RateRemap' 'Q-BLU_Day3_RateRemap'...
    'RAT1_2018-10-05_09-42-15' 'RAT4_2019-06-15_10-46-59' 'RAT4_2019-06-19_10-00-07' 'RAT5_2019-06-24_10-43-24' 'RAT5_2019-06-26_10-27-15'};

calculate_jump_distance(folders,'common')

cd(workingDir)
calculate_jump_distance(folders,'global remapped')

cd(workingDir)
extract_replay_info_batch(folders);


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
method = {'ripple 0-3', 'ripple 3-5','ripple 5-10','ripple 10 and above'};
plot_ground_truth_compare_ripple_replay_quality(folders,option,method)


%% Jump Distance

clear all
% workingDir = 'D:\ground_truth_replay_analysis\Dropbo_data8_normalised_within';
workingDir = 'P:\ground_truth_replay_analysis\Dropbo_data8';
cd(workingDir)

% workingDir = 'P:\ground_truth_replay_analysis\Dropbo_data8';
cd(workingDir)
% method = {'wcorr','linear','path','spearman'};
% method = {'spearman'};
option = 'common';
% cd('P:\ground_truth_replay_analysis\Dropbo_data8 - copy');
folders = {'2019-06-20_09-55-42' 'N-BLU_Day7_Ctrl-16x30_no_reexp' 'N-BLU_Day8_Ctrl-15-NoRest-15' 'Q-BLU_Day2_RateRemap' 'Q-BLU_Day3_RateRemap'...
    'RAT1_2018-10-05_09-42-15' 'RAT4_2019-06-15_10-46-59' 'RAT4_2019-06-19_10-00-07' 'RAT5_2019-06-24_10-43-24' 'RAT5_2019-06-26_10-27-15'};

method = {'jump distance 1','jump distance 0.6','jump distance 0.4','jump distance 0.2'};
plot_ground_truth_compare_jump_distance(folders,option,method)

% method = {'jump distance 1-0.6','jump distance 0.6-0.4','jump distance 0.4-0.2','jump distance 0.2 and below'};
% % Discarded because each event is decoded against two different tracks and
% % hence two different maximal jump distance. We 
% % plot_ground_truth_jump_distance_replay_quality(folders,option,method)

%% The effect of shuffle
clear all
% workingDir = 'D:\ground_truth_replay_analysis\Dropbo_data8_normalised_within';
workingDir = 'P:\ground_truth_replay_analysis\Dropbo_data8';
cd(workingDir)


% cd('P:\ground_truth_replay_analysis\Dropbo_data8 - copy');
folders = {'2019-06-20_09-55-42' 'N-BLU_Day7_Ctrl-16x30_no_reexp' 'N-BLU_Day8_Ctrl-15-NoRest-15' 'Q-BLU_Day2_RateRemap' 'Q-BLU_Day3_RateRemap'...
    'RAT1_2018-10-05_09-42-15' 'RAT4_2019-06-15_10-46-59' 'RAT4_2019-06-19_10-00-07' 'RAT5_2019-06-24_10-43-24' 'RAT5_2019-06-26_10-27-15'};
option = 'common';
% method = {'wcorr','spearman','linear','path'};

% method = {'time shuffle','place shuffle','place field shuffle','spike shuffle'};
method = {'spike shuffle','place field shuffle','place shuffle','time shuffle'};
plot_ground_truth_compare_single_shuffle(folders,option,method)


% method = {'three shuffles','two PRE shuffles','PRE place POST place','PRE place POST time','two POST','place bin circular shift'};
method = {'place bin circular shift','wcorr 1 shuffle + jump distance','two POST','PRE place POST time','PRE place POST place','two PRE shuffles','three shuffles'};
plot_ground_truth_compare_shuffles(folders,option,method)


%% Method optimisation

clear all
% workingDir = 'D:\ground_truth_replay_analysis\Dropbo_data8_normalised_within';
workingDir = 'P:\ground_truth_replay_analysis\Dropbo_data8';
cd(workingDir)

% folders = { '2019-06-20_09-55-42' 'N-BLU_Day7_Ctrl-16x30_no_reexp' 'N-BLU_Day8_Ctrl-15-NoRest-15' 'Q-BLU_Day2_RateRemap' 'Q-BLU_Day3_RateRemap'};
folders = {'2019-06-20_09-55-42' 'N-BLU_Day7_Ctrl-16x30_no_reexp' 'N-BLU_Day8_Ctrl-15-NoRest-15' 'Q-BLU_Day2_RateRemap' 'Q-BLU_Day3_RateRemap'...
    'RAT1_2018-10-05_09-42-15' 'RAT4_2019-06-15_10-46-59' 'RAT4_2019-06-19_10-00-07' 'RAT5_2019-06-24_10-43-24' 'RAT5_2019-06-26_10-27-15'};
option = 'common';
% method = {'wcorr','spearman','linear'};
method = {'wcorr 1 shuffle','wcorr 1 shuffle + jump distance','wcorr 2 shuffles','wcorr 3 shuffles'...
    'spearman median spike','spearman all spikes',...
    'linear 1 shuffle','linear 2 shuffles'};
plot_ground_truth_optimisation(folders,option,method)

methods = {'wcorr 1 shuffle','wcorr 1 shuffle + jump distance','wcorr 2 shuffles','wcorr 3 shuffles'...
    'spearman median spike','spearman all spikes',...
    'linear 1 shuffle','linear 2 shuffles'};
plot_replay_score_vs_log_odds_difference(folders,option,methods)
%% May 2023 plot

% folders = { '2019-06-20_09-55-42' 'N-BLU_Day7_Ctrl-16x30_no_reexp' 'N-BLU_Day8_Ctrl-15-NoRest-15' 'Q-BLU_Day2_RateRemap' 'Q-BLU_Day3_RateRemap'};
folders = {'2019-06-20_09-55-42' 'N-BLU_Day7_Ctrl-16x30_no_reexp' 'N-BLU_Day8_Ctrl-15-NoRest-15' 'Q-BLU_Day2_RateRemap' 'Q-BLU_Day3_RateRemap'...
    'RAT1_2018-10-05_09-42-15' 'RAT4_2019-06-15_10-46-59' 'RAT4_2019-06-19_10-00-07' 'RAT5_2019-06-24_10-43-24' 'RAT5_2019-06-26_10-27-15'};
option = 'common';
% method = {'wcorr','spearman','linear'};
method = {'wcorr 1 shuffle','wcorr 1 shuffle + jump distance','wcorr 2 shuffles','wcorr 3 shuffles'...
    'spearman median spike','spearman all spikes',...
    'linear 1 shuffle','linear 2 shuffles'};

plot_FDR_adjusted_log_odds_difference(folders,option,method)

%% Within normalisation vs across-track normalisation
clear all
% workingDir = 'P:\ground_truth_replay_analysis\Dropbo_data8_normalised_within';
workingDir = 'P:\ground_truth_replay_analysis\Dropbo_data8';
cd(workingDir)

% folders = { '2019-06-20_09-55-42' 'N-BLU_Day7_Ctrl-16x30_no_reexp' 'N-BLU_Day8_Ctrl-15-NoRest-15' 'Q-BLU_Day2_RateRemap' 'Q-BLU_Day3_RateRemap'};
folders = {'2019-06-20_09-55-42' 'N-BLU_Day7_Ctrl-16x30_no_reexp' 'N-BLU_Day8_Ctrl-15-NoRest-15' 'Q-BLU_Day2_RateRemap' 'Q-BLU_Day3_RateRemap'...
    'RAT1_2018-10-05_09-42-15' 'RAT4_2019-06-15_10-46-59' 'RAT4_2019-06-19_10-00-07' 'RAT5_2019-06-24_10-43-24' 'RAT5_2019-06-26_10-27-15'};
option = 'common';
method = {'wcorr normalised within','wcorr normalised across','linear normalised within','linear normalised across'};
% plot_ground_truth_compare_test(folders,option,method)
plot_ground_truth_compare_normalisation(folders,option,method)

%% log odd distribution at different p values
clear all
load log_odd_wcorr
states = [-1,0,1,2];
p_val_threshold = [-1.3 -2 -2.3];
% p_val_threshold = [-1.3 -2 -2.3 -3;
track_id = [];
track_id = log_odd.track;
epoch_track_index = [];

for p = 1:length(p_val_threshold)
    for k=1:length(states) % sort track id and behavioral states (pooled from 10 sessions)
        state_index = intersect(find(log_odd.behavioural_state==states(k)),...
            find(log_odd.pvalue <= p_val_threshold(p)));

        if k == 1 % PRE
            epoch_track_index{p}{1}{1} = intersect(state_index,find(track_id==1));
            epoch_track_index{p}{2}{1} = intersect(state_index,find(track_id==2));
        elseif k == 2 % POST
            epoch_track_index{p}{1}{3} = intersect(state_index,find(track_id==1));
            epoch_track_index{p}{2}{3} = intersect(state_index,find(track_id==2));
        elseif k == 3 % RUN Track 1
            epoch_track_index{p}{1}{2} = intersect(state_index,find(track_id==1));
        elseif k == 4 % RUN Track 2
            epoch_track_index{p}{2}{2} = intersect(state_index,find(track_id==2));
        end
    end
end


for epoch = 1:3
    epoch_index{p}{epoch} = [epoch_track_index{p}{1}{epoch} epoch_track_index{p}{2}{epoch}];
end

data = log_odd.common_zscored.original(2,:);


for p = 1:3
    fig = figure(1)
    fig.Position = [680 300 575 575];
    filename = sprintf('log odd distribution (p <= %.1i)',10^(p_val_threshold(p)));
    sgtitle(filename);
    for epoch = 1:3
        subplot(2,2,epoch)
        scatter(data(epoch_track_index{p}{1}{epoch}),0.7+rand(1,length((epoch_track_index{p}{1}{epoch})))*(1.3-0.7),3,'r','filled')
        hold on
        scatter(data(epoch_track_index{p}{2}{epoch}),1.7+rand(1,length((epoch_track_index{p}{2}{epoch})))*(2.3-1.7),3,'b','filled')
        if epoch == 1
            title('PRE')
        elseif epoch == 2
            title('RUN')

        elseif epoch == 3
            title('POST')
        end
        xlabel('z-scored log odd')
        yticks([1,2])
        yticklabels({'Track 1','Track 2'})
        ax = gca;
        set(ax,'LineWidth',1.5)
        ax.YAxis.TickDirection =  'out';       %repeat for XAxis
        ax.YAxis.TickLength =  [.005 1];       %repeat for XAxis
        ax.XAxis.TickDirection =  'out';       %repeat for XAxis
        ax.XAxis.TickLength =  [.005 1];       %repeat for XAxis
        ax.FontSize = 12;
    end

    cd Figure
    %     set(gcf,'PaperSize',[20 10]); %set the paper size to what you want
    %     print(gcf,'filename','-dpdf') % then print it
    filename = sprintf('log odd distribution (p value %.1i).pdf',10^(p_val_threshold(p)));
    saveas(gcf,filename)
    cd ..
    clf
end