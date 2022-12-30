function CSC=extract_CSC(varargin)
% INPUTS:
%   - 'all','hpc', 'ofc' or 'vc'
% LOADS:
%   - best CSC for each type of oscillation 
%   - dropped_samples
% all filtering is FIR with a set to 1 (b varies with frequency range)

parameters= list_of_parameters;
CSC_filter_length= parameters.CSC_filter_length; %default 11 with filtfilt
load('extracted_dropped_samples.mat');
    
if ~isempty(varargin)
    switch varargin{1}
        case 'all'
            MODE = [1 1];
        case 'hpc'
            MODE = [1 0];
        case 'ofc'
            MODE = [0 1];
        case 'vc'
            MODE = [0 1];
    end
else
    MODE= [1 1];
end

% Find the best CSC for each type of oscillation 
load('best_CSC')
ch = [];
ch_labels = [];
if MODE(1)==1
    ch = [ch, best_channels.bestCSC_theta  best_channels.bestCSC_ripple best_channels.bestCSC_delta best_channels.bestCSC_ThetaRippleDiff_NormalizedMethod];
    ch_labels = [ch_labels,{'best_theta' } {'best_ripple'} {'best_delta'} {'best_theta_low_ripple'}];
elseif  MODE(2)==1
    ch = [ch, best_channels.bestCSC_spindle  best_channels.bestCSC_high_gamma];
    ch_labels = [ch_labels, {'best_spindle'} {'best_high_gamma'}];
    
end 

if exist('raw')==7
    cd raw
    LFP_files=dir('*.ncs');
    dupl_LFP_files= dir('CSC*_*.ncs');
    cd ..
else
    LFP_files=dir('*.ncs');
    dupl_LFP_files= dir('CSC*_*.ncs');
end

% check if there are multiple .ncs files  
if ~isempty(dupl_LFP_files)
    CSC1_i= find(contains({LFP_files.name},'CSC1.ncs'));
    CSC_dupl_i= find(contains({dupl_LFP_files.name},'CSC1_'));
    multiples_CSC= [LFP_files(CSC1_i), dupl_LFP_files(CSC_dupl_i)];
    [~,largest_file]= max([multiples_CSC.bytes]);
    extension= strsplit(multiples_CSC(largest_file).name,'_');
    disp(['selecting _'  extension{2} ' files...']);
    LFP_files= LFP_files([LFP_files.bytes] == max(multiples_CSC.bytes));
end


c=1;
for jj=1:length(ch)  
    CSC(jj).channel=ch(jj);
    CSC(jj).channel_label=ch_labels{jj};
    if ~isempty(dupl_LFP_files)
        filename= LFP_files(contains({LFP_files.name},['CSC' num2str(ch(jj)) '_'])).name;
    else
        filename=strcat(['CSC' num2str(ch(jj)) '.ncs']);
    end
    SR = get_sample_rate(filename);
    [data,timestamps,~] = read_neuralynx_file(filename);
    disp(['extracting CSC' num2str(ch(jj))])
    Samples = data.Samples(:);  %Flatten
    timestamps = 1e-6*timestamps; %convert from microseconds to seconds
    Timestamps_buffer = linspace(0,(511/SR),512);
    Time = NaN(size(Samples));
    for i=1:512
        Time(i:512:end) = timestamps+Timestamps_buffer(i);
    end
    
    % Correct for dropped records/frames
    Time(dropped_samples) = [];
    Samples(dropped_samples) = [];
    CSCtime =(min(Time):(1/SR):max(Time))';
    signal_interp = interp1(Time,Samples,CSCtime,'linear');
    clear Time; clear Samples; clear timestamps; clear data;
    
    % parameters 
    CSC(jj).time_scale  = 'seconds';
    CSC(jj).original_SR = SR;
    CSC(jj).SR = 1000;
    CSC(jj).downsampling_factor = round(SR/CSC(jj).SR);
    CSC(jj).CSCraw  = resample(signal_interp,1,CSC(jj).downsampling_factor);
    CSC(jj).CSCtime = downsample(CSCtime,CSC(jj).downsampling_factor);
    clear signal_interp;  clear CSCtime;
    CSC(jj).filename = filename;
    
    switch cell2mat(ch_labels(jj))
        case {'best_theta','best_theta_low_ripple'}               
            
            disp('            filtering theta...');
            % theta filter         
            filter_type  = 'bandpass';
            filter_width = parameters.theta_filter;                 % range of frequencies in Hz you want to filter between
            filter_order = round(6*CSC(jj).SR/(max(filter_width)-min(filter_width)));  % creates filter with length 0.417 s for theta
            norm_freq_range = filter_width/(CSC(jj).SR/2); % SR/2 = nyquist freq i.e. highest freq that can be resolved
            b_theta = fir1(filter_order, norm_freq_range,filter_type);
            CSC(jj).theta = filtfilt(b_theta,1,CSC(jj).CSCraw);
            CSC(jj).theta_zscore = zscore(abs(hilbert(CSC(jj).theta)));
            
        case 'best_ripple'

            % ripple filter
            disp('            filtering ripples...');
            filter_type  = 'bandpass';
            filter_width = parameters.ripple_filter; %[125 300]
            filter_order = round(6*CSC(jj).SR/(max(filter_width)-min(filter_width)));    
            norm_freq_range = filter_width/(CSC(jj).SR/2); 
            b_ripple = fir1(filter_order, norm_freq_range,filter_type);
            CSC(jj).ripple = filtfilt(b_ripple,1,CSC(jj).CSCraw);
            CSC(jj).ripple_zscore = zscore(smooth(abs(hilbert(CSC(jj).ripple)),15)'); 
            
        case 'best_spindle'
            % spindle filter
            disp('            filtering spindles...');
            filter_type  = 'bandpass';
            filter_width = parameters.spindle_filter; 
            filter_order = round(6*CSC(jj).SR/(max(filter_width)-min(filter_width)));   % creates filter with length 294 s for spindles               
            norm_freq_range = filter_width/(CSC(jj).SR/2); 
            b_spindle = fir1(filter_order, norm_freq_range,filter_type);
            CSC(jj).spindle = filtfilt(b_spindle,1,CSC(jj).CSCraw);
            CSC(jj).spindle_zscore = zscore(abs(hilbert(CSC(jj).spindle)));
            
            % fast oscillation filter
            disp('            filtering fast oscillations...');
            filter_type  = 'bandpass';
            filter_width = parameters.fast_osc_filter; 
            filter_order = round(6*CSC(jj).SR/(max(filter_width)-min(filter_width)));   %           
            norm_freq_range = filter_width/(CSC(jj).SR/2); 
            b_fast = fir1(filter_order, norm_freq_range,filter_type);
            CSC(jj).fast_osc = filtfilt(b_fast,1,CSC(jj).CSCraw);
            CSC(jj).fast_osc_zscore = zscore(abs(hilbert(CSC(jj).fast_osc)));
            
        case 'best_delta'
            
            % delta filter
            disp('            filtering delta...');
            filter_type  = 'bandpass';
            filter_width = parameters.delta_filter; 
            filter_order = round(6*CSC(jj).SR/(max(filter_width)-min(filter_width)));   % creates filter with length 0.05 s for high gamma               
            norm_freq_range = filter_width/(CSC(jj).SR/2); 
            b_delta = fir1(filter_order, norm_freq_range,filter_type);
            CSC(jj).delta = filtfilt(b_delta,1,CSC(jj).CSCraw);
            CSC(jj).delta_zscore = zscore(abs(hilbert(CSC(jj).delta)));
    
        case 'best_high_gamma'
            % high gamma filter
            disp('            filtering high gamma...');
            filter_type  = 'bandpass';
            filter_width = parameters.high_gamma_filter; 
            filter_order = round(6*CSC(jj).SR/(max(filter_width)-min(filter_width)));   % creates filter with length 0.05 s for high gamma               
            norm_freq_range = filter_width/(CSC(jj).SR/2); 
            b_h_gamma = fir1(filter_order, norm_freq_range,filter_type);
            CSC(jj).high_gamma = filtfilt(b_h_gamma,1,CSC(jj).CSCraw);
            CSC(jj).high_gamma_zscore = zscore(abs(hilbert(CSC(jj).high_gamma)));
            
            % fast oscillation filter
            disp('            filtering fast oscillations...');
            filter_type  = 'bandpass';
            filter_width = parameters.fast_osc_filter; 
            filter_order = round(6*CSC(jj).SR/(max(filter_width)-min(filter_width)));   %           
            norm_freq_range = filter_width/(CSC(jj).SR/2); 
            b_fast = fir1(filter_order, norm_freq_range,filter_type);
            CSC(jj).fast_osc = filtfilt(b_fast,1,CSC(jj).CSCraw);
            CSC(jj).fast_osc_zscore=zscore(abs(hilbert(CSC(jj).fast_osc)));
    
            % low gamma filter
            disp('            filtering low gamma...');
            filter_type  = 'bandpass';
            filter_width = parameters.low_gamma_filter; 
            filter_order = round(6*CSC(jj).SR/(max(filter_width)-min(filter_width)));   % creates filter with length 0.05 s for high gamma               
            norm_freq_range = filter_width/(CSC(jj).SR/2); 
            b_l_gamma = fir1(filter_order, norm_freq_range,filter_type);
            CSC(jj).low_gamma = filtfilt(b_l_gamma,1,CSC(jj).CSCraw);
            CSC(jj).low_gamma_zscore = zscore(abs(hilbert(CSC(jj).low_gamma)));
    end
    
end

save('extracted_CSC','CSC','-v7.3');
end
