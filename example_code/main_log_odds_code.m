%% Example code for calculating z-scored log odds
% Please make sure to add example_code folder into the path or cd to the
% path

% First, we will create simulated place cells and simulated linear replay
% trajectories for two tracks
number_of_place_cells = 50;
no_of_events = 100;
[place_fields,t,replay_time,replay_spikes,ground_truth_replay]=generate_linear_replay_events(number_of_place_cells,no_of_events);

% Then, we will calculate log odds and then z-score relative to the track label
% ratemap shuffled distribution
shuffles = 1000;

track = 1;
opposite_track = 2;
bin_size = 0.02;
zscored_log_odds = [];
% This process would take 30-60 mins for 200 events
tic
parfor event = 1:10
    tic
    % Calculate Original log odds
    [~,log_odds]=calculate_log_odds(ground_truth_replay(event).spikes,place_fields,bin_size, track, opposite_track, 1);
    
    % Calculate shuffled distribution
    [~,log_odds_shuffled]=calculate_log_odds(ground_truth_replay(event).spikes,place_fields,bin_size, track, opposite_track, 1000);
    
    % Calculate z-scored log odds
    zscored_log_odds(event) = (log_odds-mean(log_odds_shuffled))/std(log_odds_shuffled);
    toc
end
zscored_log_odds(isinf(zscored_log_odds)) = nan;
save zscored_log_odds zscored_log_odds
toc

% Mean log odds based on ground truth (100% accuracy) vs mean log odds
% based on random event section
log_odds_diff = nanmean(zscored_log_odds(cat(ground_truth_replay.track_id)==1))...
    - nanmean(zscored_log_odds(cat(ground_truth_replay.track_id)==2));

s = RandStream('mrg32k3a','Seed',1); % Set random seed for resampling
random_id = randperm(s,length(ground_truth_replay));
random_id1 = random_id(1:length(random_id)/2);
random_id2 = random_id(1+length(random_id)/2:end);

log_odds_diff_random = nanmean(zscored_log_odds(random_id1))...
    - nanmean(zscored_log_odds(random_id2));

fig = figure;
fig.Position = [1100 145 560 810]
bar([log_odds_diff log_odds_diff_random])
xticklabels({'Based on ground truth','Randomly selected events'})
ylabel('Mean ')
xtickangle(45)
ax = gca;
ax.FontSize = 12;
set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
title({'Comparision between mean log odds difference' ...
    ' between Track 1 sequence and Track 2 sequence' ...
    'when selected based on ground truth or randomly selected'})

