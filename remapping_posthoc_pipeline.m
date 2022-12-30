% RATE REMAPPING _ POST HOC ANALYSIS PIPELINE

% get stable cells first/second half of laps, adds to subset_of_cells
get_stable_cells(folders,'wcorr');
rate_remapping_TRACK_PAIRS(folders,'wcorr');
plot_rate_remapping_NEW('subset','stable cells laps','epoch',{'PRE','POST'},...
                                                'x_label','peak in-field rate difference (Hz)','y_label','peak replay rate difference (Hz)');
plot_rate_remapping_NEW('subset','stable cells laps','epoch',{'sleepPRE','sleepPOST'});

% other variables
plot_rate_remapping_NEW('subset','stable cells laps','epoch',{'PRE','POST'},...
                                                'y_var',{'replay_spike_diff_nonZero'},...
                                                'x_label','peak in-field rate difference (Hz)','y_label','mean number of spikes per event');
plot_rate_remapping_NEW('subset','stable cells laps','epoch',{'PRE','POST'},...
                                                'y_var',{'median_replay_rate_diff'},...
                                                'x_label','peak in-field rate difference (Hz)','y_label','median replay rate (Hz)');
plot_rate_remapping_NEW('subset','stable cells laps','epoch',{'PRE','POST'},...
                                                'x_var',{'plfield_AUC_diff'},...
                                                'x_label','place field AUC difference','y_label','peak replay rate difference (Hz)');
% single peak cells
plot_rate_remapping_NEW('subset',{'stable cells laps','single place field cells'},'epoch',{'PRE','POST'},...
                                                'x_label','peak in-field rate difference (Hz)','y_label','peak replay rate difference (Hz)');
% rate modulated cells
plot_rate_remapping_NEW('subset',{'stable cells laps','rate modulated cells'},'epoch',{'PRE','POST'},...
                                                'x_label','peak in-field rate difference (Hz)','y_label','peak replay rate difference (Hz)');
% Controls
within_field_analysis(method) % finds place cells with multiple fields, adds output to Remapping main structure
remapping_classification_bootstrap(folders) % finds rate modulated cells
plot_control_track_peakFR() % increase track PeakFR threshold
plot_correlation_sleep_evolution(folders,'wcorr', 'POST','events')
plot_rate_remapping_NEW('subset','stable cells laps','epoch',{'awakePOST','sleepPOST'},...
                                                'x_label','peak in-field rate difference (Hz)','y_label','peak replay rate difference (Hz)');
plot_correlation_RUN_evolution(folders, 'wcorr', 'cumulative_laps')
                                           
% Plotting
plot_rate_remapping_NEW(varargin) %plots correlations
get_distributions_tracks_replay(varagin) % NEEDS UPDATE - plots FR distributions for tracks & replay (?)
summary_data(folders,option) % summary table
num_events = calculate_number_replay_events;

get_corr_tracks(folders,'method','pop_vector');

% GLMs/Regressions
prepare_data_for_regression(varargin); % under construction
plot_rate_remapping_NEW('subset','predictive cells','epoch',{'PRE','POST'})

% Examine REPLAY firing distributions for each cell in analysis, for significant track1,2 and
% all replay events
get_firing_stats_replay(folders,method);
compare_track_FR % Quick comparison of peak FR distributions between tracks

% Awake replay on track (RUN)
plot_correlation_RUN_SLEEP('epochs',{'PRE','POST'},'target_epoch',{'RUN'},'vars','mean_max_FR_replay_diff'); %main
plot_correlation_RUN_SLEEP('epochs',{'PRE','POST'},'target_epoch',{'RUN'},'vars','mean_max_FR_replay_diff','scoring','spearman'); %spearman
plot_correlation_RUN_SLEEP('epochs',{'PRE','POST'},'target_epoch',{'RUN'},'vars','mean_max_FR_replay_diff','control','rate_detection_control'); % track peak rate shuffle
plot_correlation_RUN_SLEEP('epochs',{'PRE','POST'},'target_epoch',{'RUN'},'vars','mean_max_FR_replay_diff','control','replay_rate_shuffle_detection'); % replay rate shuffle

% get examples of raw data with couple cells highlighted
example_neuron_raw_data('cell_id',[5097 5029]);
plot_confusion_matrices;

%%% SHUFFLES %%%
plot_intrinsic_bias_shuffle_distribution(shuffle_type)

% SPEARMAN
rate_remapping_TRACK_PAIRS(folders,'spearman');
plot_rate_remapping_NEW('x_var',{'place_field_diff'},'y_var',{'mean_max_FR_replay_diff'},'epochs',{'PRE','POST'},...
                                                            'scoring','spearman','subset','stable cells laps',...
                                                            'x_label','peak in-field rate difference (Hz)','y_label','peak replay rate difference (Hz)');

% TRACK PEAK RATE SHUFFLE
track_peak_rate_shuffle(folders,'rate_intrinsic_bias');
track_peak_rate_shuffle(folders,'rate_detection_control');
plot_intrinsic_bias_shuffle_distribution('rate_bias');
plot_rate_remapping_NEW('x_var',{'place_field_diff'},'y_var',{'mean_max_FR_replay_diff'},'epochs',{'PRE','POST'},...
                                                            'subset','stable cells laps','control','rate_detection_control',...
                                                            'x_label','peak in-field rate difference (Hz)','y_label','peak replay rate difference (Hz)');

% REPLAY RATE SHUFFLE
replay_rate_shuffle_analysis(folders,'replay_rate_shuffle');
replay_rate_shuffle_analysis(folders,'replay_rate_shuffles_only');
plot_intrinsic_bias_shuffle_distribution('replay_rate_bias');
plot_rate_remapping_NEW('x_var',{'place_field_diff'},'y_var',{'mean_max_FR_replay_diff'},'epochs',{'PRE','POST'},...
                                                            'subset','stable cells laps','control','replay_rate_shuffle_detection',...
                                                            'x_label','peak in-field rate difference (Hz)','y_label','peak replay rate difference (Hz)');

                                                        
%% individual sessions
for i=1:length(folders)
rat_folder= folders(i,1);
i
remapping= rate_remapping_TRACK_PAIRS(rat_folder,'wcorr',0);
plot_rate_remapping_NEW('use_mat',remapping,'subset','stable cells laps',...
                        'epoch',{'PRE','POST'},...
                        'x_label','peak in-field rate difference (Hz)',...
                        'y_label','peak replay rate difference (Hz)');
% ax= get(gcf,'children');
% set([ax(:)],'plotboxaspectratio',[1.3 1 1]);
% ylim([ax(:)],[-3 3]); yticks([ax(:)],[-3 0 3]); yticklabels([ax(:)],[-3 0 3]);
% xlim([ax(:)],[-30 30]); xticks([ax(:)],[-30:15:30]); xticklabels([ax(:)],[-30:15:30]);
% saveas(gcf,['..\Figures\Single session correlations\session' num2str(i) '.pdf']);
close(gcf);
end