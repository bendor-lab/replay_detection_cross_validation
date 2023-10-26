# replay_detection_cross_validation
Matlab analysis code for validating hippocampal replay detection methods. For more information about the replay analysis framework, please read: https://www.biorxiv.org/content/10.1101/2022.12.12.520040v3

## Prerequisites and installation
- MATLAB 2019a or newer (working on MATLAB 2021a and MATLAB 2022b).
- MATLAB Toolboxes: Parallel Computing Toolbox Signal Processing Toolbox Statistics and Machine Learning Toolbox
- Add this repository into the path by using 'set path' -> 'Add with subfolders' (installation time less than 10 seconds)

## Code for main analysis


## Code for eLife revision
Two main codes are uploaded for all the additional analyses in response to eLife revision feedback
- **elife_replay_model_BH_test.mat** for testing BH method using simple replay model
- **elife_add_on_analyses.mat** for the rest of the additional analyses including
  - Creating two additional randomized dataset (cross experiment cell id randomization and place field circular shifted randomization)
  - Extracting log odds and replay information
  - Validating the use of cell-id randomized dataset
  - Plotting replay detection performance within each session
  - Plotting the effect of ripple threshold on replay detection
  - Plotting 4 individual examples of replay events

 ## Sample Code for Bayesian decoding and log odds calculation
 Go to **main_log_odds_code.mat** to run an example code for Bayesian decoding and log odds quantification based on Poisson-spiking place cell and simulated linear replay trajectories
- 
