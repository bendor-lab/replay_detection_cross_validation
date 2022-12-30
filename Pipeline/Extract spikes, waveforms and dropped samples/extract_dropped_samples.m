function extract_dropped_samples(varargin)
%extracts timestamps of samples dropped from Neuralynx file during the recording
%note that these samples must be removed at a later step in the pipeline,
%as during "drops" samples are drawn from the existing samples in the buffer and do not
%reflect real data.

if isempty(varargin)
    d= dir('raw/*.ncs');
    d1_i= find(contains({d.name},'CSC1.ncs'));
    d2= dir('raw/CSC1_*.ncs');
    if ~isempty(d2)
        d3= [d(d1_i), d2];
        % check which file is the largest
        disp('multiple .ncs files available, selecting largest...');
        [~,largest_file]= max([d3.bytes]);
        filename= d3(largest_file).name;
        disp(['using ' filename]);
    else
        disp('automatically selecting CSC1.ncs');
        filename= 'CSC1.ncs';
    end
else
    filename= varargin{1};
    disp(['Using ' filename]);
end

SR=get_sample_rate(filename);
[data,timestamps,~]=read_neuralynx_file(filename);

Samples=data.Samples(:);  %Flatten
timestamps=1e-6*timestamps; %convert from microseconds to seconds
Timestamps_buffer=linspace(0,(511/SR),512);
Time=NaN(size(Samples));
for i=1:512
    Time(i:512:end)=timestamps+Timestamps_buffer(i);
end
dropped_frames=find(data.NumberOfValidSamples<512);
dropped_samples=[];
for i=1:length(dropped_frames)
    dropped_samples=[dropped_samples ((dropped_frames(i)-1)*512+((1+data.NumberOfValidSamples(dropped_frames(i))):512))];
end
error_times=Time(dropped_samples);

save('extracted_dropped_samples.mat','dropped_samples','error_times','Time','-v7.3');
end