function all_PSD=extract_PSD(varargin)
% calculates PSD for selected channels
% INPUT:
%   - vector of channels
%   - no input: all channels will have their PSD calculated
%   - use sleep states or not
%   e.g: extract_PSD(channels,sleep) or extract_PSD(channels,[])
% LOADS:
%   - position data
%   - sleep state (based on velocity thresholding)

load('extracted_position.mat');
%constants for PSD
new_SR   = 1000;            %new sampling rate
nfft = 1024;            % number of point in FT
win  = hanning(nfft);   % hanning window
noverlap = [];          % amount of overlap of sections of signal ([] = default = 50% overlap)

if length(varargin)>=2 && strcmp(varargin{2},'sleep')==1
    % take advantage of sleep states to calculate different PSDs awake
    [states,state_name]=get_awake_and_sleep_times;
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

for c = 1:length(LFP_files)
    LFP_filename = LFP_files(c).name;
    if isempty(dupl_LFP_files)
        name_split = strsplit(LFP_filename,{'.','CSC'},'CollapseDelimiters',true);
        CSC_number = str2num(name_split{2});
        ncs_name= [{'CSC'} {'.ncs'}];
    elseif ~isempty(dupl_LFP_files)
        name_split = strsplit(LFP_filename,{'_','CSC'});
        CSC_number = str2num(name_split{2});
        ncs_name= [{'CSC'} {['_' name_split{3}]}];
    else
        disp('something went wrong in the filenames');
        keyboard;
    end

    available_channels(c) = CSC_number;
end

available_channels=unique(available_channels);
if ~isempty(varargin{1})
    available_channels= intersect(varargin{1},available_channels);
end
available_channels(find(isnan(available_channels)))=[];  %remove channels that are NaNs
if isempty(available_channels)
    disp('ERROR- no available channels to perform PSD')
end
for q = 1:length(available_channels)
    %get filename and Sample Rate
    filename = strcat([ncs_name{1} num2str(available_channels(q)) ncs_name{2}])
    original_SR = get_sample_rate(filename);
    if ~isempty(original_SR)
        %to save time do not reinterpolate signal
        all_PSD(q).original_SR = original_SR;
        all_PSD(q).SR = new_SR;
        all_PSD(q).downsampling_factor = round(original_SR/all_PSD(q).SR);
        if length(varargin)>=2 && strcmp(varargin{2},'sleep')==1
            for state_id = 1:length(states)
                state = states{state_id};
                for f = 1:size(state,1)
                    [data,~,~] = read_neuralynx_file(filename, [1 1 1 1 1],1,4,1e6*position.t(state(f,1:2)));  %extract CSC files
                    Samples = data.Samples(:);  %Flatten
                    CSC = resample(Samples,1,all_PSD(q).downsampling_factor); %resample to 1000 Hz
                    [Pxx,F] = pwelch(CSC, win, noverlap, nfft, new_SR);    %calculate PSD, constants at beginning of function
                    all_PSD(q).(['PSD_F_' state_name{state_id}]) = F;
                    all_PSD(q).(['PSD_' state_name{state_id}])(f,:) = 10*log10(Pxx); % to convert to dB/Hz if required
                    all_PSD(q).CSCchannel = available_channels(q);
                    all_PSD(q).filename=filename;
                end
            end
        else
            [data,~,~] = read_neuralynx_file(filename);  %extract CSC files
            Samples = data.Samples(:);  %Flatten
            steps = floor(length(Samples)/(original_SR*60*15)); %15 min blocks
            index = 1:(original_SR*60*15);
            for f=1:(steps-1)
                CSC=resample(Samples(index+((f-1)*max(index))),1,all_PSD(q).downsampling_factor);
                [Pxx,F] = pwelch(CSC, win, noverlap, nfft, new_SR);    %calculate PSD, constants at beginning of function
                all_PSD(q).PSD_F= F;
                all_PSD(q).PSD(f,:)=log10(Pxx); % to convert to dB/Hz if required
            end
            all_PSD(q).CSCchannel=available_channels(q);
            all_PSD(q).filename=filename;
        end
    end
end

PSD_states.states = states;
PSD_states.states_name = state_name;
save('PSD_data','all_PSD','PSD_states');
end
