function parameters=list_of_parameters

%for place field calculations
% parameters.speed_threshold=.5;  %should be 2 cm/s
% parameters.place_field_smoothing=21; 
% parameters.min_smooth_peak=1;
% parameters.min_raw_peak=1;
% parameters.max_mean_rate=1;

parameters.MUA_filter_length= 41; % in samples (ms)
parameters.MUA_filter_alpha= 2;

% for LFP analysis
parameters.theta_filter= [4 12]; % in Hz
parameters.ripple_filter= [125 300];
parameters.spindle_filter= [9 17];%check for HVS (7-12) - other papers 7-14/12-15/10-15
parameters.delta_filter= [1 4]; %other papers - 0.5-2 slow delta & 2-4 fast delta
parameters.fast_osc_filter= [100 300];
parameters.high_gamma_filter= [40 100]; %other papers 35-100
parameters.low_gamma_filter= [17 40];

% for tracking
parameters.max_pixel_distance = 50;
parameters.max_pixel_distance_jump = 40; % before 100
parameters.position_filter_length = 101; %51
parameters.speed_filter_length = 25; %1s for filtfilt
parameters.max_distance_jump = 40;
parameters.max_speed_sleepbox = 100;

% For place field calculation
parameters.speed_threshold = 4;  %should be 2 cm/s
parameters.speed_threshold_laps = 5; 
parameters.speed_threshold_max = 50; 
parameters.place_field_smoothing = 10; % multiply by x_bins_width to have cm
parameters.place_field_smoothing_bayesian = 2; % multiply by x_bins_width to have cm
parameters.min_smooth_peak = 0.5;
parameters.min_raw_peak = 1;
parameters.max_mean_rate = 5;   
parameters.x_bins_width = 2;
parameters.x_bins_width_bayesian = 10;
parameters.half_width_threshold = 0.0005; %50us

% For waveform extraction
parameters.nSamplesForSpike = [-24, 24];
parameters.SR = 30000; %30kHz
parameters.nChannels = 4;

% For CSC extraction
parameters.CSC_filter_length = 15; %default 11 with filtfilt

%for replay detection
parameters.min_zscore=0;
parameters.max_zscore=3;
parameters.max_search_length=300;
parameters.min_event_duration=0.1; %80ms
parameters.max_event_duration=1; %1s

% For bayesian decoding
parameters.replay_bin_width = 0.02;
parameters.run_bin_width = 0.25;
parameters.position_bin_width = 10;
parameters.smoothing_number_of_bins = 5;
parameters.bayesian_threshold = 1e-4;

parameters.colors.purple = [0.49,0.18,0.56];
parameters.colors.orange = [0.93,0.69,0.13];
parameters.colors.brick = [0.7647,0.0784,0.1137];
parameters.colors.blue = [0, 0.4470, 0.7410];
parameters.plot_color_line= {'r','b','g','c','m','y','k'};
parameters.plot_color_dashed_line= {'r--','b--','g--','c--','m--','y--','k--'};
parameters.plot_color_symbol= {'ro','bx','g+','cs','md','y*','k^'};
parameters.plot_color_dot= {'r.','b.','g.','c.','m.','y.','k.'};  
parameters.legend={'track 1','track 2','track 3','track 4','track 5','track 6'};

parameters.min_react_duration = 0.05;
parameters.max_event_duration = 0.750;
parameters.min_event_duration = 0.1;
parameters.max_event_duration = 0.750;


parameters.rng_seed_remapping= [11 12 13 14 15];
end 