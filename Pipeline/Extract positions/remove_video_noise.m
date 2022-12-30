function positions=remove_video_noise(targets, timestamps)
% DB_Nov '18
% Uses targets provided by Nlyx to plot all the x and y detected by the
% software during recording (including both LEDs light position and other
% sources of noise picked up by the software).


positions.x=[];
positions.y=[];
positions.timestamps=[];

for i=1:1000:length(targets)
    j=i+999;
    if j>length(targets)
        j=length(targets);
    end

    [x,y,color,valid_targets] = ExtractFromTargets(targets(:,i:j));
    t=timestamps(i:j);

   for k=1:max(valid_targets)
        index=find(valid_targets>=k);
        positions.x=[positions.x; x(index,k)];
        positions.y=[positions.y; y(index,k)];
        positions.timestamps=[positions.timestamps; t(index)'];
   end

end