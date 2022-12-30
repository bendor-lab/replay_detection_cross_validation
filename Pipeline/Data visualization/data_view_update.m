function data_view_update(input)
global time
global play
global data
global t
global FIG
global GUI_handle
global replay_index;
global parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(input,'load') %load data
    parameters=list_of_parameters;
    replay_index=[];
    GUI_handle=gcbf;
    FIG=figure;
    if exist('extracted_dropped_samples.mat')==2
        disp('loading dropped samples')
        load('extracted_dropped_samples.mat');
        data.error_times=error_times;
    else
        data.error_times=[];
    end
    if exist('extracted_position.mat')==2
        disp('loading position file')
        load('extracted_position.mat');
        data.position=position;
        time.rezero_time=extract_injection_time;
        if isempty(time.rezero_time)
            time.rezero_time=min(position.t);
        end
        data.position.t=data.position.t-time.rezero_time;
        data.position.sleepbox=data.position.sleepbox-time.rezero_time;
        data.position.x=rescale(data.position.x);
        data.position.y=rescale(data.position.y);
        data.position.v=clip(data.position.v,0,5)/5;
        for i=1:length(data.position.linear)
            data.position.linear(i).linear=data.position.linear(i).length*rescale(data.position.linear(i).linear);
        end
    else
        data.position=[];
    end

    if exist('extracted_sleep_state.mat')==2
        load('extracted_sleep_state.mat');
        data.sleep_state=sleep_state;
        data.sleep_state.time=data.sleep_state.time-time.rezero_time;
    else
        data.sleep_state=[];
    end
    
    if exist('MUA_clusters.mat')==2
        disp('loading multiunit_activity file')
        load('MUA_clusters.mat');
        data.mua= MUA;
        data.mua.time_bins=data.mua.time_bins-time.rezero_time;
        data.mua.zscore=clip(data.mua.zscored,0,10)/10;
%         load('multiunit_activity.mat');
%         data.multiunit_activity=multiunit_activity;
%         data.multiunit_activity.timestamp=data.multiunit_activity.timestamp-time.rezero_time;
    else
        data.mua=[];
    end
    if exist('extracted_clusters.mat')==2
        disp('loading clusters file')
        load('extracted_clusters.mat');
        data.clusters=clusters;
        data.clusters.spike_times=data.clusters.spike_times-time.rezero_time;

    else data.clusters=[];
    end
    if exist('extracted_place_fields.mat')==2
        disp('loading place_fields')
        load('extracted_place_fields.mat');
        data.place_fields=place_fields;
    else data.clusters=[];
    end
    
    
    if exist('extracted_CSC.mat')==2
        disp('loading CSC file')
        load('extracted_CSC.mat');
        % this is if you have multiple brain areas
        CSC_idx= find(~cellfun(@isempty,{CSC.ripple}));
        data.CSC= CSC(CSC_idx);
        data.CSC.CSCtime=data.CSC.CSCtime-time.rezero_time;
        data.CSC.CSCraw= 4*rescale(data.CSC.CSCraw)-2;
        data.CSC.ripple_zscore= clip(data.CSC.ripple_zscore,0,10)/10;
    else data.CSC=[];   
    end
        if exist('extracted_replay_events.mat')==2
        disp('loading replay events')
        load('extracted_replay_events.mat');
        data.replay=replay;
        data.replay.onset=data.replay.onset-time.rezero_time;
        data.replay.offset=data.replay.offset-time.rezero_time; 
        data.replay_events=make_dash(data.replay.onset,data.replay.offset);
        
%         data.replay_events_mua=make_dash(replay_mua.onset,replay_mua.offset);
%         data.replay_events_mua=data.replay_events_mua-time.rezero_time;
%         data.replay_events_ripple=make_dash(replay_ripple.onset,replay_ripple.offset);
%          data.replay_events_ripple=data.replay_events_ripple-time.rezero_time;
        
        else data.replay=[];
        end
    if exist('estimated_position.mat')==2
        disp('loading decoded position')
        load('estimated_position.mat');
        %load('estimated_position_concat.mat')
        data.estimated_position=estimated_position;
       else
        data.estimated_position=[];
    end
    
    time.min_time=min(data.position.t);
    time.max_time=max(data.position.t);
    time.start_experiment_time=time.min_time;
    time.stop_experiment_time=time.max_time;
    time.all_time=time.max_time-time.min_time;
    write_to_GUI('min_time',num2str(time.min_time));
    write_to_GUI('max_time',num2str(time.max_time));
    write_to_GUI('all_time',num2str(time.all_time));
    time.time_scale=time.all_time;
elseif strcmp(input,'zoom') %ZOOM
    time.zoom_index=read_value_from_GUI('zoom_index');
    time.all_time=time.time_scale*time.zoom_index;
    time.max_time=time.min_time+time.all_time;
    write_to_GUI('min_time',num2str(time.min_time));
    write_to_GUI('max_time',num2str(time.max_time));
    write_to_GUI('all_time',num2str(time.all_time));
elseif strcmp(input,'scroll') %SCROLL
    time.scroll_index=read_value_from_GUI('scroll_index');
    time.min_time=time.start_experiment_time+time.scroll_index*time.time_scale;
    time.max_time=time.min_time+time.all_time;
    write_to_GUI('min_time',num2str(time.min_time));
    write_to_GUI('max_time',num2str(time.max_time));
    write_to_GUI('all_time',num2str(time.all_time));
elseif strcmp(input,'new_figure');
    FIG=figure;
else
    check_times;
    if strcmp(input,'plot')
        time.min_time=read_number_from_GUI('min_time');
        time.max_time=read_number_from_GUI('max_time');
        time.all_time=read_number_from_GUI('all_time');
        plot_data(2);
        replay_index=[];
    elseif strcmp(input,'next_replay')
        time.min_time=read_number_from_GUI('min_time');
        if isempty(replay_index)
            replay_index=find((data.replay.onset-time.min_time)>0);
            replay_index=replay_index(1);
            time.min_time=pseudo_round(data.replay.onset(replay_index)-1);
        else
            replay_index=replay_index+1;
            if replay_index<length(data.replay.onset)
                time.min_time=pseudo_round(data.replay.onset(replay_index)-1);
            end
        end
          write_to_GUI('min_time',num2str(time.min_time));
          write_to_GUI('max_time',num2str(time.min_time+2));
          write_to_GUI('all_time',num2str(2));
           time.min_time=read_number_from_GUI('min_time');
        time.max_time=read_number_from_GUI('max_time');
        time.all_time=read_number_from_GUI('all_time');
        plot_data(1);
    elseif strcmp(input,'play')
        replay_index=[];
        play=read_value_from_GUI('play');
         time.all_time=read_number_from_GUI('all_time');
        if(play)
             write_to_GUI('play','PLAY');
             if isempty(timerfind)
                 t = timer;
                 t.TimerFcn= @plot_timer;
                 t.BusyMode='drop';
                 if time.all_time<300
                     t.Period = .25;
                 elseif time.all_time<1000
                    t.Period = 0.5; 
                 else
                    t.Period = 1;  
                 end
                 t.ExecutionMode = 'fixedRate';
                 start(t)
             end
        else
            write_to_GUI('play','STOP');
            if ~isempty(timerfind)
                stop(t)
                delete(t)
            end
        end
    end
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=pseudo_round(x)
y=round(10*x)/10;
end

function out=read_value_from_GUI(tag)
global GUI_handle
EditHandle = findobj(GUI_handle,'Tag',tag);
out=get(EditHandle,'Value');
end

function out=read_number_from_GUI(tag)
global GUI_handle
EditHandle = findobj(GUI_handle,'Tag',tag);
out=str2num(get(EditHandle,'String'));
end

function out=read_file_from_GUI(tag)
global GUI_handle
EditHandle = findobj(GUI_handle,'Tag',tag);
out=get(EditHandle,'String');
end

function write_to_GUI(tag,data)
global GUI_handle
f1=findobj(GUI_handle,'Tag',tag);
set(f1,'string',data);
end

function enable_GUI(tag,choice)
global GUI_handle
f1=findobj(GUI_handle,'Tag',tag);
if choice==1
    set(f1,'enable','on');
elseif choice==-1
    set(f1,'enable','off');
end
end

function set_GUI(tag,value)
global GUI_handle
f1=findobj(GUI_handle,'Tag',tag);
set(f1,'value',value);
end

function visible_GUI(tag,choice)
global GUI_handle
f1=findobj(GUI_handle,'Tag',tag);
if choice==1
    set(f1,'visible','on');
elseif choice==-1
    set(f1,'visible','off');
end
end

function check_times
global time;
time.min_time=pseudo_round(read_number_from_GUI('min_time'));
time.max_time=pseudo_round(read_number_from_GUI('max_time'));
time.all_time=pseudo_round(read_number_from_GUI('all_time'));
if time.all_time<0.1
    time.all_time=0.1;
end
%update scroll and zoom
if time.min_time<time.start_experiment_time
    time.min_time=time.start_experiment_time;
end
time.max_time=time.min_time+time.all_time;
if time.max_time>time.stop_experiment_time
    time.max_time=time.stop_experiment_time;
end
time.zoom_index=time.all_time/time.time_scale;
time.scroll_index=(time.min_time-time.start_experiment_time)/time.time_scale;
if time.zoom_index>1
    time.zoom_index=1;
end
if time.zoom_index<0
    time.zoom_index=0.001;
end
if time.scroll_index>1
    time.scroll_index=1;
end
if time.scroll_index<0
    time.scroll_index=0;
end
set_GUI('scroll_index',time.scroll_index);
set_GUI('zoom_index',time.zoom_index);
write_to_GUI('min_time',num2str(pseudo_round(time.min_time)));
write_to_GUI('max_time',num2str(pseudo_round(time.max_time)));
write_to_GUI('all_time',num2str(pseudo_round(time.all_time)));
end

function x=clip(x,a,b)
x(find(x<a))=a;
x(find(x>b))=b;
end

function y=make_dash(a,b)
 y=NaN(1,3*length(a));
        y(1:3:end)=a;
        y(2:3:end)=b;
end

function raster_plot(x,y,c,h)
x2(1:3:length(x)*3)=x;
x2(2:3:length(x)*3)=x;
x2(3:3:length(x)*3)=NaN;
y2(1:3:length(x)*3)=y;
y2(2:3:length(x)*3)=y+h;
y2(3:3:length(x)*3)=NaN;
if isempty(c)
    plot(x2,y2);
else
    plot(x2,y2,'Color',c,'LineWidth',2);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






function plot_data(option)
global data
global FIG
global time;
global parameters
    %PLOT POSITION

    figure(FIG);
     ax(1)=subplot('Position',[.1 0.05 0.85 .60]);
    hold off
    plot([time.min_time time.max_time],[0 0],'w.')
    refresh;
    hold on
   xlabel('time(s)');
    y.bottom=0.5;  y.top=1.5; y.text{1}=' ';
      p_index=find(data.position.t>time.min_time & data.position.t<time.max_time);
    
      
      if option~=2
   y.text{1}='sleepbox';
     index=find(data.position.sleepbox>time.min_time & data.position.sleepbox<time.max_time);
    if ~isempty(index)
     plot(data.position.sleepbox(index),ones(size(index)),'m.')
    end
    
    index=find(data.error_times>time.min_time & data.error_times<time.max_time);
    if ~isempty(index)
     plot(data.error_times(index),ones(size(index)),'kx')
    end
    if ~isempty(data.sleep_state)
          index=find(data.sleep_state.time>time.min_time & data.sleep_state.time<time.max_time);
        plot(data.sleep_state.time(index),0.5+data.sleep_state.state(index),'c')
    end
   

   
   y.bottom(end+1)=y.top(end)+0.5;
   y.top(end+1)=y.bottom(end)+1;
   y.text{end+1}='speed';
   
  
    plot(data.position.t(p_index),y.bottom(end)+data.position.v(p_index),'k')

 
   y.bottom(end+1)=y.top(end)+0.5;
   y.top(end+1)=y.bottom(end)+1;
   y.text{end+1}='x/y position';
   
    plot(data.position.t(p_index),y.bottom(end)+data.position.x(p_index),'r') 
    plot(data.position.t(p_index),y.bottom(end)+data.position.y(p_index),'g')
   


    y.bottom(end+1)=y.top(end)+1.5;
   y.top(end+1)=y.bottom(end)+1;
   y.text{end+1}='MUA';
    
     %plot mua
     index=find(data.mua.time_bins>time.min_time & data.mua.time_bins<time.max_time); 
      plot(data.mua.time_bins(index),y.bottom(end)+data.mua.zscore(index),'k');
     
  index=find(data.CSC.CSCtime>time.min_time & data.CSC.CSCtime<time.max_time); 
  plot(data.CSC.CSCtime(index),y.bottom(end)+-data.CSC.ripple_zscore(index),'r');
    
  
  index=find(data.replay_events>time.min_time & data.replay_events<time.max_time);
  plot(data.replay_events(min(index):max(index)),y.bottom(end)+ones(size(min(index):max(index))),'b','LineWidth',4)

%    index=find(data.replay_events_mua>time.min_time & data.replay_events_mua<time.max_time);
%   plot(data.replay_events_mua(min(index):max(index)),.2+y.bottom(end)+ones(size(min(index):max(index))),'k','LineWidth',2)
%  
%    index=find(data.replay_events_ripple>time.min_time & data.replay_events_ripple<time.max_time);
%   plot(data.replay_events_ripple(min(index):max(index)),.1+y.bottom(end)+ones(size(min(index):max(index))),'r','LineWidth',2)
 end
 

 
          %plot csc
        y.bottom(end+1)=y.top(end)+1.5;
   y.top(end+1)=y.bottom(end)+1;
 y.text{end+1}='LFP';
  index=find(data.CSC.CSCtime>time.min_time & data.CSC.CSCtime<time.max_time); 
      plot(data.CSC.CSCtime(index),y.bottom(end)+data.CSC.CSCraw(index),'k');
      data.index=index;
  

 c=parameters.plot_color_line;
     spike_index=find(data.clusters.spike_times>time.min_time & data.clusters.spike_times<time.max_time);
     spike_times=data.clusters.spike_times(spike_index)';
     spike_id=data.clusters.spike_id(spike_index);
    for i=1:length(data.position.linear)
        y.bottom(end+1)=y.top(end)+0.5;
        number_of_units=max(data.place_fields.track(i).sorted_good_cells);
        y.top(end+1)=y.bottom(end)+8;
        y.text{end+1}=strcat(['track ' num2str(i) ' cells']);
        index=[];
        spike_id_sorted=[];
        spike_id_unique=[];
        spike_times_sorted=[];
         spike_times_unique=[];
        for j=1:length(data.place_fields.track(i).sorted_good_cells)
            index=find(spike_id==data.place_fields.track(i).sorted_good_cells(j));
            spike_id_sorted=[spike_id_sorted j*ones(1,length(index))];
            spike_times_sorted=[spike_times_sorted spike_times(index)];
          
        end
        raster_plot(spike_times_sorted,y.bottom(end)+10*spike_id_sorted/number_of_units,c{i},10/number_of_units);
    end
%  
if option~=2
y.bottom(end+1)=y.top(end)+0.5;
        number_of_units=max(data.place_fields.other_cells);
        y.top(end+1)=y.bottom(end)+6;
        y.text{end+1}='other cells';
        index=[];
        spike_id_sorted=[];
        spike_times_sorted=[];
        [temp,sort_index]=sort(data.place_fields.track(1).mean_rate_track(data.place_fields.other_cells));
        for j=1:length(data.place_fields.other_cells)
            index=find(spike_id==data.place_fields.other_cells(sort_index(j)));
            spike_id_sorted=[spike_id_sorted j*ones(1,length(index))];
            spike_times_sorted=[spike_times_sorted spike_times(index)];
        end
        raster_plot(spike_times_sorted,y.bottom(end)+10*spike_id_sorted/number_of_units,'k',10/number_of_units);
end

    y.bottom(end+1)=y.top(end)+.5;
   y.top(end+1)=y.bottom(end)+.5;
     y.text{end+1}='';
    %plot clusters
     for i=1:length(y.top)
        plot([time.min_time time.max_time],[y.top(i)+0.25 y.top(i)+0.25],'k--');
    end
    axis([time.min_time time.max_time y.bottom(1) y.top(end)]);
    yticklabels(y.text)
    yticks((y.bottom+y.top)/2)
   
    %%%%%%%%%%%%plotting bayesian decoding
reverse_bone=flipud(colormap(bone));
number_of_tracks=length(data.position.linear);
plot_height=.3/number_of_tracks;

for i=1:number_of_tracks
    ax(1+i)=subplot('Position',[0.1 .65+(i-1)*plot_height 0.85 plot_height]);
    hold off
   
    if time.all_time<=10 & ~isempty(data.estimated_position)
        index1=find(data.estimated_position(i).replay_time_centered>(time.min_time +time.rezero_time) & data.estimated_position(i).replay_time_centered<(time.max_time+time.rezero_time));
        estimated_position= data.estimated_position(i).replay(:,index1);
      lower_color_limit=(1./size(estimated_position,1));
        imagesc(data.estimated_position(i).replay_time_centered(index1)-time.rezero_time,data.estimated_position(i).position_bin_centres, estimated_position,[lower_color_limit 0.5]);
        hold on
        colormap(reverse_bone);
         plot3(data.position.t(p_index),100*data.position.linear(i).linear(p_index),ones(1,length(p_index)),'Color',c{i},'linewidth',2)
    elseif time.all_time<=300 & ~isempty(data.estimated_position)
        index1=find(data.estimated_position(i).run_time_centered>(time.min_time +time.rezero_time) & data.estimated_position(i).run_time_centered<(time.max_time+time.rezero_time));
        estimated_position= data.estimated_position(i).run(:,index1);
          lower_color_limit=(1./size(estimated_position,1));
        imagesc(data.estimated_position(i).run_time_centered(index1)-time.rezero_time,data.estimated_position(i).position_bin_centres, estimated_position,[lower_color_limit 0.5]);
        hold on
        colormap(reverse_bone);
        plot3(data.position.t(p_index),100*data.position.linear(i).linear(p_index),ones(1,length(p_index)),'Color',c{i},'linewidth',2)
    else
        plot(data.position.t(p_index),100*data.position.linear(i).linear(p_index),'Color',c{i},'linewidth',2);
        hold on
        plot([time.min_time time.max_time],[0 0],'w.')
    end
    axis([time.min_time time.max_time min(data.estimated_position(i).position_bin_centres)-10 max(data.estimated_position(i).position_bin_centres)+10]);
    axis xy;
    ylabel(strcat(['track' num2str(i)]));
end
refresh;
linkaxes(ax,'x');

end



function plot_timer(mTimer,~)
global t
global time
%update scroll
time.scroll_index=read_value_from_GUI('scroll_index');
play_speed=read_number_from_GUI('play_speed');
time.min_time=time.min_time+play_speed*time.all_time/10;
time.max_time=time.min_time+time.all_time;
time.scroll_index=(time.min_time-time.start_experiment_time)/time.time_scale;

time.min_time=pseudo_round(time.start_experiment_time+time.scroll_index*time.time_scale);
time.max_time=pseudo_round(time.min_time+time.all_time);
time.all_time=pseudo_round(time.all_time);
write_to_GUI('min_time',num2str(time.min_time));
write_to_GUI('max_time',num2str(time.max_time));
write_to_GUI('all_time',num2str(time.all_time));
plot_data(0);
if (time.scroll_index>=1 | time.max_time>=time.stop_experiment_time)
    time.scroll_index=1;
    set_GUI('play',0);
    write_to_GUI('play','STOP');
    stop(t);
    delete(t);
end
set_GUI('scroll_index',time.scroll_index);

end
