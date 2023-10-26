# replay_detection_cross_validation
Matlab analysis code for validating hippocampal replay detection methods. For more information about the replay analysis framework, please read: https://www.biorxiv.org/content/10.1101/2022.12.12.520040v3

## Prerequisites and installation
- MATLAB 2019a or newer (working on MATLAB 2021a and MATLAB 2022b).
- MATLAB Toolboxes: Parallel Computing Toolbox Signal Processing Toolbox Statistics and Machine Learning Toolbox
- Add this repository into the path by using 'set path' -> 'Add with subfolders' (installation time less than 10 seconds)

## Preprocessing 
Go to **batch_analysis.m** in Pipeline folder for most of the basic preprocessing pipelines including:

1 - spikes & video ;
2 - position & laps ;
3 - sleep processing ;
4 - CSC & PSD ;
5 - Place cells and Bayesian decoding ;
6 - Candidate replay event extraction and replay decoding ;
7 - replay scoring and significance

## Main analysis
**main_analysis_code.mat** is used for the majority of the analysis 
- **ground_truth_replay_sequence_analysis.mat** for sequence-based replay detection and quantification
- **ground_truth_log_odd_analysis.mat** for creating cell-id randomized dataset as well as sequenceless decoding (log odds quantification) for both original and randomized dataset
- **calculate_jump_distance.mat** for calculating jump distance
- **extract_replay_info_batch.mat** for extracting key replay information necessary for subsequent analysis
- **plot_ground_truth_compare_ripple_replay_quality.mat** for analysing the effect of ripple power on replay detection
- **plot_ground_truth_compare_jump_distance.mat** for analysing the effect of jump distance on replay detection
- **plot_ground_truth_compare_single_shuffle.mat** for analysing the effect of different shuffling proceure on replay detection
- **plot_ground_truth_compare_shuffles.mat** for analysing the effect of adding stricter criteria (i.e. jump distance or more shuffling procedures) on replay detection
- **plot_ground_truth_optimisation.mat** for analysing the effect of different detection methods (including weighted correlation, linear fitting and spearman) on replay detection
- **plot_ground_truth_compare_normalisation.mat** for analysing the effect of cross-track and within track normalisation on replay detection

## Code for eLife revision
Two main codes are uploaded for all the additional analyses in response to eLife revision feedback
- **elife_replay_model_BH_test.mat** for testing BH method using a simple replay model
- **elife_add_on_analyses.mat** for the rest of the additional analyses including
  - Creating two additional randomised dataset (cross experiment cell-id randomisation and place field circular shifted randomisation)
  - Extracting log odds and replay information
  - Validating the use of cell-id randomized dataset
  - Plotting replay detection performance within each session
  - Plotting the effect of ripple threshold on replay detection
  - Plotting 4 individual examples of replay events

 ## Sample Code for Bayesian decoding and log odds calculation
 Go to **main_log_odds_code.mat** to run an example code for Bayesian decoding and log odds quantification based on Poisson-spiking place cell and simulated linear replay trajectories.
- All dependent functions are in **example_code** folder:
  - **generate_linear_replay_events.mat** for creating simulated place cells with spike trains that follow the Poisson process and for generating simulated linear replay trajectories for both tracks.
  - **calculate_log_odds.mat** for running Bayesian decoding and log odds quantification
-Follow log odds quantification, the original log odds will be z-scored relative to the shuffled log odds distribution (track label ratemap shuffle, see paper for more information).
- As a showcase of the utility of mean log odds difference, we compare the mean log odds difference when track 1 and track 2 replay events are based on ground truth (true track identity) and are randomly selected. You will see the mean log odds difference between track 1 and track 2 sequences when selected correctly is substantially larger than when randomly selected.
- With modification, anyone can use this code to analyse real replay data.
