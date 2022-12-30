function see_replay_update(input)

global time
global replay_index
global data
global FIG
global GUI_handle
global parameters
global tagged_replay
global new_index
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(input,'load') %load data
    parameters = list_of_parameters;
    GUI_handle = gcbf;
    FIG = figure;
    tagged_replay = [];
    new_index  = [];
    data=[];
    %%% Extract Position data %%%
    if exist('extracted_position.mat')==2
        disp('loading position file')
        load('extracted_position.mat');
        data.position = position;
        
        time.rezero_time = min(data.position.t);
        
        data.position.t = data.position.t-time.rezero_time;
        
    else
        data.position = [];
    end
    
    %%% Extract MUA data %%%
   
    if exist('MUA_clusters.mat')==2  % from updated MUA code
        disp('loading MUA clusters')
        load('MUA_clusters.mat');
        
        data.mua.time_bins = MUA.time_bins-time.rezero_time;
        data.mua.zscore = rescale(MUA.zscored);
    else
        data.mua = [];
    end
    
    %%% Extract Place fields info %%%
    load extracted_place_fields_BAYESIAN;
    data.place_fields_BAYESIAN=place_fields_BAYESIAN;
    disp('loading place fields')
    
    %%% Extract cluster information (spikes and cluster ID) %%%
    if exist('extracted_clusters.mat')==2
        disp('loading clusters file')
        load('extracted_clusters.mat');
        data.clusters = clusters;
        data.clusters.spike_times = data.clusters.spike_times-time.rezero_time;
        clear clusters;     
    else
        data.clusters = [];
    end
    
    %%% Extract CSC information %%%
    if exist('extracted_CSC.mat')==2
        disp('loading CSC file')
        load('extracted_CSC.mat');        
        % this is if you have multiple brain areas
        hpc_idx = find(~cellfun(@isempty,{CSC.ripple}));
        theta_idx = find(~cellfun(@isempty,{CSC.theta}));
        data.CSC.CSCtime = CSC(hpc_idx).CSCtime - time.rezero_time;
        data.CSC.CSCraw = rescale(CSC(hpc_idx).CSCraw)+0.05;
        data.CSC.ripple_zscore = rescale(CSC(hpc_idx).ripple_zscore'); %scale from 0 to 1
        data.CSC.theta_zscore = rescale(CSC(theta_idx(1)).theta_zscore);
    else
        data.CSC=[];
    end
    
    %%% Extract replay event features %%%
    if exist('extracted_replay_events.mat')==2
        disp('loading replay events')
        load('extracted_replay_events.mat');
        data.replay = replay;
        data.replay.onset = data.replay.onset-time.rezero_time;
        data.replay.offset = data.replay.offset-time.rezero_time;
        data.replay.thresh_adapted_end = replay.thresh_adapted_end;
        data.replay.thresh_adapted_start = replay.thresh_adapted_start;
        data.replay.midpoint=replay.midpoint-time.rezero_time;
        replay_index = 1;        
    else
        data.replay = [];
    end
    
    if exist('scored_replay.mat')==2
        disp('loading scored replay events')
        load('scored_replay.mat');
        for track=1:length(scored_replay)
            for event=1:length(scored_replay(track).replay_events)
                data.replay_score.linear_score(track,event)=scored_replay(track).replay_events(event).linear_score;
                data.replay_score.weighted_corr_score(track,event)=scored_replay(track).replay_events(event).weighted_corr_score;
                data.replay_score.spearman_score(track,event)=scored_replay(track).replay_events(event).spearman_score;
                data.replay_score.path_score(track,event)=scored_replay(track).replay_events(event).path_score;
            end
        end
        if exist('significant_replay_events.mat')==2
            load('significant_replay_events.mat');
            data.significant_replay_events = significant_replay_events;
            data.significant_replay_events.score=zeros(length(significant_replay_events.track),length(significant_replay_events.pre_ripple_threshold_index));
            data.significant_replay_events.event_segment_best_score=zeros(length(significant_replay_events.track),length(significant_replay_events.pre_ripple_threshold_index));
              data.significant_replay_events.bayesian_bias=zeros(length(significant_replay_events.track),length(significant_replay_events.pre_ripple_threshold_index));
        for track=1:length(significant_replay_events.track)
                for event=1:length(significant_replay_events.track(track).replay_score)
                     data.significant_replay_events.score(track,significant_replay_events.track(track).index(event))=significant_replay_events.track(track).replay_score(event);
                     data.significant_replay_events.event_segment_best_score(track,significant_replay_events.track(track).index(event))=data.significant_replay_events.track(track).event_segment_best_score(event);
                     data.significant_replay_events.bayesian_bias(track,significant_replay_events.track(track).index(event))=data.significant_replay_events.track(track).bayesian_bias(event);
               end
            end
        else
            data.significant_replay_events=[];
            significant_replay_events.pre_ripple_threshold_index = 1:length(data.replay.onset);
        end
    else
        data.replay_score=[];
    end
    
    data.replay.onset = data.replay.onset(significant_replay_events.pre_ripple_threshold_index);
    data.replay.offset = data.replay.offset(significant_replay_events.pre_ripple_threshold_index);
    data.replay.thresh_adapted_end = data.replay.thresh_adapted_end(significant_replay_events.pre_ripple_threshold_index);
    data.replay.thresh_adapted_start = data.replay.thresh_adapted_start(significant_replay_events.pre_ripple_threshold_index);
    data.replay.midpoint = data.replay.midpoint(significant_replay_events.pre_ripple_threshold_index);
    
    data.replay.zscore = data.replay.zscore(significant_replay_events.pre_ripple_threshold_index);
    data.replay.ripple_peak = data.replay.ripple_peak(significant_replay_events.pre_ripple_threshold_index);
    data.replay.duration = data.replay.duration(significant_replay_events.pre_ripple_threshold_index);
    data.replay.speed = data.replay.speed(significant_replay_events.pre_ripple_threshold_index);
    data.replay.spike_count = data.replay.spike_count(significant_replay_events.pre_ripple_threshold_index);
    data.replay.neuron_count = data.replay.neuron_count(significant_replay_events.pre_ripple_threshold_index);
    data.replay.mean_theta = data.replay.mean_theta(significant_replay_events.pre_ripple_threshold_index);   
  
    data.replay_score.spearman_score= data.replay_score.spearman_score(significant_replay_events.pre_ripple_threshold_index); 
    data.replay_score.linear_score=data.replay_score.linear_score(significant_replay_events.pre_ripple_threshold_index); 
    data.replay_score.weighted_corr_score=data.replay_score.weighted_corr_score(significant_replay_events.pre_ripple_threshold_index); 
    data.replay_score.path_score=data.replay_score.path_score(significant_replay_events.pre_ripple_threshold_index); 
    
    
     
    if exist('decoded_replay_events.mat')==2
        disp('loading decoded replay events')
        load('decoded_replay_events.mat');
        for track=1:length(decoded_replay_events)
            for event=1:length(decoded_replay_events(track).replay_events(significant_replay_events.pre_ripple_threshold_index))
                data.decoded_position{track,event}=modify_decoding(decoded_replay_events(track).replay_events(significant_replay_events.pre_ripple_threshold_index(event))); %.decoded_position;         
            end
        end
    end
    
   
    % Set time 
    time.min_time = min(data.position.t);
    time.max_time = max(data.position.t);
    time.start_experiment_time = time.min_time;
    time.stop_experiment_time  = time.max_time;
    time.all_time = time.max_time-time.min_time;
    time.time_scale = time.all_time;
    
    % Sort replay events by different features and save sorted indices
     new_index(1,:) = 1:length(data.replay.zscore); %unsorted
    [~,new_index(2,:)] = sort(data.replay.zscore); 
    [~,new_index(3,:)] = sort(data.replay.ripple_peak);
    [~,new_index(4,:)] = sort(data.replay.duration);
    [~,new_index(5,:)] = sort(data.replay.speed);
    [~,new_index(6,:)] = sort(data.replay.spike_count);
    [~,new_index(7,:)] = sort(data.replay.neuron_count);
    [~,new_index(8,:)] = sort(data.replay.mean_theta);   
    [~,new_index(9,:)] = sort(max(data.replay_score.spearman_score,[],1)); 
    [~,new_index(10,:)] = sort(max(data.replay_score.linear_score,[],1));
    [~,new_index(11,:)] = sort(max(data.replay_score.weighted_corr_score,[],1));
    [~,new_index(12,:)] = sort(max(data.replay_score.path_score,[],1));    
      [~,new_index(13,:)] = sort(max(data.significant_replay_events.score,[],1));
   
    
    
elseif strcmp(input,'previous')  % to show previous 6 replay events
    replay_index = replay_index-6;
    if replay_index<1
        replay_index = 1;
    end
    set_GUI('slider1',replay_index/length(data.replay.onset));
    write_to_GUI('replay_count',strcat([num2str(replay_index) ' / ' num2str(length(data.replay.onset))]));
    plot_data(replay_index);

elseif strcmp(input,'next')
    replay_index = replay_index+6;  % to show next 6 replay events
    if replay_index>(length(data.replay.onset)-6)
        replay_index=(length(data.replay.onset)-6);
    end
    set_GUI('slider1',replay_index/length(data.replay.onset));
    write_to_GUI('replay_count',strcat([num2str(replay_index) ' / ' num2str(length(data.replay.onset))]));
    plot_data(replay_index);

elseif strcmp(input,'slider')
    replay_index = round(read_value_from_GUI('slider1')*length(data.replay.onset));
    if replay_index<1
        replay_index=1;
    end
    if replay_index>(length(data.replay.onset)-6)
        replay_index=(length(data.replay.onset)-6);
    end
    write_to_GUI('replay_count',strcat([num2str(replay_index) ' / ' num2str(length(data.replay.onset))]));
    plot_data(replay_index);

elseif strcmp(input,'analyze')  % to plot all good and bad replay events according to different features
    figure
    subplot(3,3,1)
    plot(data.replay.speed,data.replay.zscore,'k.')
    hold on
    plot(data.replay.speed(tagged_replay),data.replay.zscore(tagged_replay),'ro')
    xlabel('speed')
    ylabel('mua zscore')
    
    subplot(3,3,2)
    plot(data.replay.speed,data.replay.ripple_peak,'k.')
    hold on
    plot(data.replay.speed(tagged_replay),data.replay.ripple_peak(tagged_replay),'ro')
    xlabel('speed')
    ylabel('ripple_peak')
    
    subplot(3,3,3)
    plot(data.replay.speed,data.replay.duration,'k.')
    hold on
    plot(data.replay.speed(tagged_replay),data.replay.duration(tagged_replay),'ro')
    xlabel('speed')
    ylabel('duration')
    
    subplot(3,3,4)
    plot(data.replay.speed,data.replay.spike_count,'k.')
    hold on
    plot(data.replay.speed(tagged_replay),data.replay.spike_count(tagged_replay),'ro')
    xlabel('speed')
    ylabel('spike_count')
    
    subplot(3,3,5)
    plot(data.replay.speed,data.replay.neuron_count,'k.')
    hold on
    plot(data.replay.speed(tagged_replay),data.replay.neuron_count(tagged_replay),'ro')
    xlabel('speed')
    ylabel('neuron count')
    
    subplot(3,3,6)
    plot(data.replay.ripple_peak,data.replay.zscore,'k.')
    hold on
    plot(data.replay.ripple_peak(tagged_replay),data.replay.zscore(tagged_replay),'ro')
    xlabel('ripple_peak')
    ylabel('mua zscore')
    
    subplot(3,3,7)
    plot(data.replay.neuron_count,data.replay.zscore,'k.')
    hold on
    plot(data.replay.neuron_count(tagged_replay),data.replay.zscore(tagged_replay),'ro')
    xlabel('neuron count')
    ylabel('mua zscore')
    
    subplot(3,3,8)
    plot(data.replay.spike_count,data.replay.zscore,'k.')
    hold on
    plot(data.replay.spike_count(tagged_replay),data.replay.zscore(tagged_replay),'ro')
    xlabel('spike count')
    ylabel('mua zscore')
    
    subplot(3,3,9)
    plot(data.replay.spike_count,data.replay.neuron_count,'k.')
    hold on
    plot(data.replay.spike_count(tagged_replay),data.replay.neuron_count(tagged_replay),'ro')
    xlabel('spike_count')
    ylabel('neuron_count')
    
elseif strcmp(input,'1')  % classify chosen replay event as bad
    tagged_replay=[tagged_replay new_index(read_value_from_GUI('analysis_choice'),replay_index)];
elseif strcmp(input,'2')
    tagged_replay=[tagged_replay new_index(read_value_from_GUI('analysis_choice'),replay_index+1)];
elseif strcmp(input,'3')
    tagged_replay=[tagged_replay new_index(read_value_from_GUI('analysis_choice'),replay_index+2)];
elseif strcmp(input,'4')
    tagged_replay=[tagged_replay new_index(read_value_from_GUI('analysis_choice'),replay_index+3)];
elseif strcmp(input,'5')
    tagged_replay=[tagged_replay new_index(read_value_from_GUI('analysis_choice'),replay_index+4)];
elseif strcmp(input,'6')
    tagged_replay=[tagged_replay new_index(read_value_from_GUI('analysis_choice'),replay_index+5)];
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

function x=clip(x,a,b) % currently not being used
x(find(x<a))=a;
x(find(x>b))=b;
end

function y=make_dash(a,b) % currently not being used
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
    plot(x2,y2,'LineWidth',1);
else
    plot(x2,y2,c,'LineWidth',2);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_data(replay_index)
global data
global FIG
global time;
global parameters
global new_index

if read_value_from_GUI('analysis_choice')==2
    write_to_GUI('score',num2str(data.replay.zscore(new_index(read_value_from_GUI('analysis_choice'),replay_index))));

elseif read_value_from_GUI('analysis_choice')==3
    write_to_GUI('score',num2str(data.replay.ripple_peak(new_index(read_value_from_GUI('analysis_choice'),replay_index))));
    
elseif read_value_from_GUI('analysis_choice')==4
    write_to_GUI('score',num2str(data.replay.duration(new_index(read_value_from_GUI('analysis_choice'),replay_index))));
    
elseif read_value_from_GUI('analysis_choice')==5
    write_to_GUI('score',num2str(data.replay.speed(new_index(read_value_from_GUI('analysis_choice'),replay_index))));
    
elseif read_value_from_GUI('analysis_choice')==6
    write_to_GUI('score',num2str(data.replay.spike_count(new_index(read_value_from_GUI('analysis_choice'),replay_index))));
    
elseif read_value_from_GUI('analysis_choice')==7
    write_to_GUI('score',num2str(data.replay.neuron_count(new_index(read_value_from_GUI('analysis_choice'),replay_index))));

elseif read_value_from_GUI('analysis_choice')==8
    write_to_GUI('score',num2str(data.replay.mean_theta(new_index(read_value_from_GUI('analysis_choice'),replay_index))));
elseif read_value_from_GUI('analysis_choice')==9
    write_to_GUI('score', num2str(data.replay_score.spearman_score(:, new_index(read_value_from_GUI('analysis_choice'),replay_index))));
elseif read_value_from_GUI('analysis_choice')==10
    write_to_GUI('score',num2str(data.replay_score.linear_score(:, new_index(read_value_from_GUI('analysis_choice'),replay_index))));

elseif read_value_from_GUI('analysis_choice')==11
    write_to_GUI('score',num2str(data.replay_score.weighted_corr_score(:, new_index(read_value_from_GUI('analysis_choice'),replay_index))));
elseif read_value_from_GUI('analysis_choice')==12
    write_to_GUI('score',num2str(data.replay_score.path_score(:, new_index(read_value_from_GUI('analysis_choice'),replay_index))));
elseif read_value_from_GUI('analysis_choice')==13
    write_to_GUI('score',num2str(data.significant_replay_events.score(:, new_index(read_value_from_GUI('analysis_choice'),replay_index))));

else
    write_to_GUI('score',num2str(replay_index));
    
end
number_of_tracks=size(data.significant_replay_events.score,1);
rows_for_plots=number_of_tracks+4;
rr = 0;
for r = new_index(read_value_from_GUI('analysis_choice'),replay_index:(replay_index+5))
    rr = rr+1;
    time.min_time = data.replay.onset(r)-0.1;
    time.max_time = data.replay.offset(r)+0.1;
    
    figure(FIG);
    
    for track=1:number_of_tracks
    %plot replay of tracks
    subplot(rows_for_plots,6,rr+(number_of_tracks-track)*6)
    imagesc(data.decoded_position{track,r},[0 1])
    axis xy
    if data.significant_replay_events.score(track, r)>0        
    label=title(strcat(['score= ' num2str(data.significant_replay_events.score(track, r))]));
     label.Color = 'red';
     if data.significant_replay_events.event_segment_best_score(track,r)==1
           label=xlabel(strcat(['seg=entire , bias=' num2str(data.significant_replay_events.bayesian_bias(track,r))]));
 
     elseif data.significant_replay_events.event_segment_best_score(track,r)==2
           label=xlabel(strcat(['seg=1st half , bias=' num2str(data.significant_replay_events.bayesian_bias(track,r))]));
 
     elseif data.significant_replay_events.event_segment_best_score(track,r)==3
         label=xlabel(strcat(['seg=2nd half , bias=' num2str(data.significant_replay_events.bayesian_bias(track,r))]));
 
     end
     label.Color = 'red';
    end
    
    if ~isempty(find(data.significant_replay_events.multi_tracks_index==r))
    label=ylabel('multiple tracks significant');
 label.Color = 'green';
    end 
    end
    colormap(flipud(colormap('bone')))
   
   
         %%%%%%%%%%%
    subplot(rows_for_plots,6,rr+(number_of_tracks+1)*6)    
    hold off;
    c = parameters.plot_color_line;
    spike_index = find(data.clusters.spike_times>time.min_time & data.clusters.spike_times<time.max_time);
    spike_times = data.clusters.spike_times(spike_index)';
    spike_id = data.clusters.spike_id(spike_index);   

    [~,~,spike_id2]=unique(spike_id);
    raster_plot(spike_times,spike_id2,'k',1); 
    hold on;
    number_of_units = max(spike_id2);
    spike_index = find(spike_times>=data.replay.onset(r) & spike_times<=data.replay.offset(r));
    spike_times = spike_times(spike_index)';
    spike_id = spike_id(spike_index);
    spike_id2 = spike_id2(spike_index);
    raster_plot(spike_times,spike_id2,'r',1);
  axis([time.min_time time.max_time 0 number_of_units]);
    yl= ylim;
    rectangle('Position',[data.replay.onset(r) yl(1) data.replay.offset(r)-data.replay.onset(r) yl(2)-yl(1)],'FaceColor',[0.5 0.5 0.5 0.5],'EdgeColor','none');
    yticks([0:10:number_of_units])
    xlabel('time (s)');
    if mod(rr,6)==1
        ylabel('neurons');
    end

    for j=1:length(data.place_fields_BAYESIAN.track)
    spike_times_sorted{j}=spike_times;
    spike_id_sorted{j}=NaN(size(spike_id));
    for i=1:length(data.place_fields_BAYESIAN.track(j).sorted_good_cells)
        a=find(spike_id==data.place_fields_BAYESIAN.track(j).sorted_good_cells(i));
        if ~isempty(a)
        spike_id_sorted{j}(a)=i;
        end
    end
    index=find(isnan(spike_id_sorted{j}));
    spike_id_sorted{j}(index)=[];
    spike_times_sorted{j}(index)=[];
    [~,~,spike_id2_sorted{j}]=unique(spike_id_sorted{j});
    
    %%%%%%%%%%%%%%%%
     subplot(rows_for_plots,6,rr+(number_of_tracks)*6) 
    hold on;
    c = ['b','m','g','r'];
     raster_plot(spike_times_sorted{j},(j-1)*number_of_units+spike_id2_sorted{j},c(j),1);
   end
    
     
    %plot clusters (spike raster plot)
    axis([time.min_time time.max_time 0 j*number_of_units]);
    yl= ylim;
    rectangle('Position',[data.replay.onset(r) yl(1) data.replay.offset(r)-data.replay.onset(r) yl(2)-yl(1)],'FaceColor',[0.5 0.5 0.5 0.5],'EdgeColor','none');
      plot([data.replay.midpoint(r) data.replay.midpoint(r)],[0 j*number_of_units],'k--');
    yticks([0:10:j*number_of_units])
    xlabel('time (s)');
    if mod(rr,6)==1
        ylabel('neurons');
    end
%%%%%%%
 
    subplot(rows_for_plots,6,rr+(number_of_tracks+2)*6)    
    hold off
    
    %plot mua
    index = find(data.mua.time_bins>time.min_time & data.mua.time_bins<time.max_time);
    plot(data.mua.time_bins(index),data.mua.zscore(index),'k');
    hold on
    xlabel('time(s)');
    index1 = find(data.CSC.CSCtime>time.min_time & data.CSC.CSCtime<time.max_time);
    plot(data.CSC.CSCtime(index1),-2+data.CSC.ripple_zscore(index1),'k');   %plot ripple zscored
       
    %plot raw csc    
    index2 = find(data.CSC.CSCtime>time.min_time & data.CSC.CSCtime<time.max_time);   
    plot(data.CSC.CSCtime(index2),5*data.CSC.CSCraw(index2),'k');
    data.index = index2;
    
    yticks([-2 0 2]);
    if mod(rr+6,12)==1
        yticklabels({'ripple','mua','raw'});
    else
        yticklabels({});
    end
    ylim([-4 4]); yl= ylim;
    axis([time.min_time time.max_time ylim]);
    rectangle('Position',[data.replay.onset(r) yl(1) data.replay.offset(r)-data.replay.onset(r) yl(2)-yl(1)],'FaceColor',[0.5 0.5 0.5 0.5],'EdgeColor','none');
    xlabel('time (s)');
    if data.replay.thresh_adapted_start(r) ~= 0
        plot(data.replay.onset(r),5*data.CSC.CSCraw(find(min(abs(data.CSC.CSCtime-data.replay.onset(r))))),'rx');
        text(data.replay.onset(r),5*data.CSC.CSCraw(find(min(abs(data.CSC.CSCtime-data.replay.onset(r)))))+1,num2str(data.replay.thresh_adapted_start(r)));
    end
    if data.replay.thresh_adapted_end(r) ~= 0
        plot(data.replay.offset(r),5*data.CSC.CSCraw(find(min(abs(data.CSC.CSCtime-data.replay.offset(r))))),'rx');
        text(data.replay.offset(r),5*data.CSC.CSCraw(find(min(abs(data.CSC.CSCtime-data.replay.offset(r)))))+1,num2str(data.replay.thresh_adapted_end(r)));
    end
         plot([data.replay.midpoint(r) data.replay.midpoint(r)],ylim,'r--');

         
    % plot speed, if rat is in pot and theta
%%%%%%%%%%%%%%%
   subplot(rows_for_plots,6,rr+(number_of_tracks+3)*6)    
    index1 = find(data.position.t>time.min_time & data.position.t<time.max_time);
    plot(data.position.t(index1),data.position.v_cm(index1),'k'); hold on; %plot speed
    plot(xlim,[5 5],'r');
    index2 = find(data.position.sleepbox>time.min_time & data.position.sleepbox<time.max_time); %plot sleepbox position
    plot(data.position.sleepbox(index2), ones(size(data.position.sleepbox(index2))),'b.'); hold on;
    index3 = find(data.replay.onset>time.min_time & data.replay.offset<time.max_time);
        plot(xlim,[data.replay.mean_theta(index3(1)) data.replay.mean_theta(index3(1))],'m--');
    
    if ~isempty(index2) && mod(rr+12,6)==1
        legend('speed','speed limit','sleep','mean theta zscore','Location','northwest');
        legend boxoff;
        ylabel('speed (cm/s)');
    elseif mod(rr+12,6)==1
        legend('speed','speed limit','mean theta zscore','Location','northwest');
        legend boxoff;
        ylabel('speed (cm/s)');
    end
    ylim([0 10]);
    xlim([time.min_time time.max_time]);
    xlabel('time (s)');
     hold off;
    
    
end
end



function modified_decoded_event = modify_decoding(events)
modified_decoded_event = events.decoded_position;
for i = 1 : length(events.timebins_edges)-1
    spikes = find(events.spikes(:,2)>=events.timebins_edges(i) & events.spikes(:,2)<events.timebins_edges(i+1),1);  %>= and < used, to match histcount function
    if isempty(spikes)
        modified_decoded_event(:,i) = zeros(size(modified_decoded_event(:,i)));
    end
end
end

