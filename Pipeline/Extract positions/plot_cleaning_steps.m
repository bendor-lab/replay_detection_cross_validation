function plot_cleaning_steps

load('extracted_video_data.mat');
load('position_data');
load('extracted_position.mat');
parameters= list_of_parameters;

number_of_tracks= length(position.linear);
%plot pre and post processing x and y positions
figure;
subplot(2,1,1);
plot(positions.x,positions.y,'Color',[0.5 0.5 0.5]); hold on;
plot(positions.x,positions.y,'r.','MarkerSize',1);
title('raw X and Y positions')
subplot(2,1,2);
plot(position.clean.x, position.clean.y,'b.');
title('X and Y after pre-processing');

%plot post processing x and y positions on linear tracks
figure;
for i= 1:number_of_tracks
    subplot(number_of_tracks,2,2*i-1);
    plot(position.clean.x(position.linear(i).cropped_indices),position.clean.y(position.linear(i).cropped_indices),'k.');
    title(['Track ' num2str(i)]);
end

%plot linear tracks with colour map representing mapping of linear position
for i= 1:number_of_tracks
    subplot(number_of_tracks,2,2*i);
    imagesc(flipud(rot90(position.linear(i).grid.index))); 
    colormap(jet)
    hold on
    plot(position.clean.x(position.linear(i).cropped_indices),position.clean.y(position.linear(i).cropped_indices),'k.')
    set(gca,'YDir','normal')
    xlim([0 max(position.x)]);
    ylim([0 max(position.y)]);
    title(strcat(num2str(parameters.max_pixel_distance),' pixels'))
end

%plot raw x position, clearn x position, and x position on track and in sleep box 
figure;
ax(1)= subplot(4,1,1);
plot(positions.timestamps,positions.x);
ax(2)=subplot(4,1,2);
plot(position.t,position.x);
ax(3)= subplot(4,1,3);
plot(position.clean.t, position.clean.x,'r'); hold on;
plot(position.sleep.t, position.sleep.x,'b');
legend({'all','sleep'});
%plot linear position on track and in sleep box 
ax(4)=subplot(4,1,4);
for i=1:length(position.linear)
    plot(position.t,position.linear(i).linear); hold on;
end
plot(position.sleepbox,ones(size(position.sleepbox)),'m.');
linkaxes(ax,'x')

end