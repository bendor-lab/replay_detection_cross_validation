function [cropped_indices,cropped_logical] = crop_in_time(data_matrix,timestamps,track_num)
% data_matrix should be arranged with each column is a variable
% corresponding timestamps

figure;
max_trace= 0;
for i=1:size(data_matrix,2)
    plot(timestamps,data_matrix(:,i)+ max_trace); hold on;
    max_trace= max(data_matrix(:,i))+ 5;
end
title([{'Select start and end points for'}; {['TRACK ' num2str(track_num)]}])

 for hh=1:2
    current_point(hh,:)= ginput(1);
    plot([current_point(hh,1) current_point(hh,1)],ylim,'r');  
 end

first_timestamp= current_point(1,1);
last_timestamp= current_point(2,1);

% cropped_timestamps= timestamps(timestamps >= first_timestamp & timestamps<=last_timestamp);
cropped_indices= find(timestamps >= first_timestamp & timestamps<=last_timestamp);
cropped_logical= (timestamps >= first_timestamp & timestamps<=last_timestamp);

end