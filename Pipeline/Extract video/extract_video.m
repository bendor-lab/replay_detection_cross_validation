function positions=extract_video(varargin)

% If varargin = 'targets', run code to find all targets (all pixels rather
% than just the position based on the largest region)
[data,timestamps,~]=read_neuralynx_file('VT1.nvt');
timestamps=1e-6*timestamps; %convert to seconds
if strcmp(varargin,'targets')
    disp('extracting targets')
    positions=remove_video_noise(data.Targets, timestamps);
    positions.targets=1;
    positions=targets_to_position(positions);  %convert targets to unique timestamps
else
    positions.timestamps = timestamps;
    positions.x = data.ExtractedX;
    positions.y = data.ExtractedY;
    positions.targets=0;
end

save extracted_video_data positions;
end