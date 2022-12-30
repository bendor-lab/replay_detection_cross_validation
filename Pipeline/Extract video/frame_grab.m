function frame_grab

if exist('raw/VT1.mpg')==2
    v = VideoReader('raw/VT1.mpg');
    extracted_video_frame = readFrame(v);
    save extracted_video_frame extracted_video_frame
else
    disp('raw/VT1.mpg does not exist');
end