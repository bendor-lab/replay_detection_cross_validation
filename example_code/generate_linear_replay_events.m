function [place_fields,t,replay_time,replay_spikes,ground_truth_replay]=generate_linear_replay_events(number_of_place_cells,no_of_events)
% This code is used to generate place fields and then simulate linear replay trajectory and generate place 
% cell spiking (Possion process) 
% 
% Masahiro Takigawa 2023


% First create place cells for two tracks
% number_of_place_cells = 100;
place_cell_centres=(0:number_of_place_cells-1)/(number_of_place_cells-1);
parameters.rate=[10 50];   %[10 100]; %[fixed rate additive random rate]
parameters.spacing=1; %0 is random, %1 is equally spaced
parameters.bandwidth=1;  %0 is random, %1 is equally spaced
parameters.speed_threshold=0.1;
% parameters.bayesian_threshold=10.^(log2(number_of_place_cells)-log2(400)); % small value multiplied to all values to get rid of zeros

place_fields.track = [];
place_fields.track(1).sorted_good_cells = 1:1:number_of_place_cells; 

for j=1:number_of_place_cells
    place_fields.track(1).x =  0:0.05:1; % Normalised distance
    place_fields.track(1).centre(j)=place_cell_centres(j);% tuned to a different location
    s = RandStream('mrg32k3a','Seed',j); % Set random seed for resampling
    place_fields.track(1).peak_rate(j)=parameters.rate(1)+parameters.rate(2)*rand(s,1);
    place_fields.track(1).rate{j}=zeros(size(place_fields.track(1).x));
    
    index=interp1(place_fields.track(1).x,1:length(place_fields.track(1).x),place_fields.track(1).centre(j),'nearest');
    bandwidth=round(length(place_fields.track(1).x)/5);
    
    place_fields.track(1).index{j}=(index-bandwidth):(index+bandwidth);
    place_fields.track(1).shape{j}=place_fields.track(1).peak_rate(j)*gausswin(length(place_fields.track(1).index{j}),5).^2;
    
    m=find(place_fields.track(1).index{j}>0 & place_fields.track(1).index{j}<=length(place_fields.track(1).x));

    if ~isempty(m)
        place_fields.track(1).rate{j}(place_fields.track(1).index{j}(m))=place_fields.track(1).shape{j}(m);
    end
end

% Create Track 2 place cell map by cell id randomisation
shuffled_id = randperm(number_of_place_cells);
place_fields.track(2) = place_fields.track(1);
place_fields.track(2).sorted_good_cells = shuffled_id;

for j=1:number_of_place_cells
    place_fields.track(2).x =  0:0.001:1;
    place_fields.track(2).centre(shuffled_id(j))=place_cell_centres(shuffled_id(j));% tuned to a different location
    place_fields.track(2).peak_rate(shuffled_id(j))=parameters.rate(1)+parameters.rate(2)*rand(1);
    place_fields.track(2).rate{shuffled_id(j)}=zeros(size(place_fields.track(2).x));
    
    index=interp1(place_fields.track(2).x,1:length(place_fields.track(2).x),place_fields.track(2).centre(j),'nearest');
    if parameters.bandwidth == 0
        bandwidth=round(length(place_fields.track(2).x)/(4+2*rand(1)));
    else 
        bandwidth=round(length(place_fields.track(2).x)/5);
    end
    
    place_fields.track(2).index{shuffled_id(j)}=(index-bandwidth):(index+bandwidth);
    place_fields.track(2).shape{shuffled_id(j)}=place_fields.track(2).peak_rate(shuffled_id(j))*gausswin(length(place_fields.track(2).index{shuffled_id(j)}),5).^2;
    
    m=find(place_fields.track(2).index{shuffled_id(j)}>0 & place_fields.track(2).index{shuffled_id(j)}<=length(place_fields.track(2).x));

    if ~isempty(m)
        place_fields.track(2).rate{shuffled_id(j)}(place_fields.track(2).index{shuffled_id(j)}(m))=place_fields.track(2).shape{shuffled_id(j)}(m);
    end
    
end

fig = figure;
fig.Position = [680 215 950 760];
c=1;
for kk=1:length(place_fields.track)
    for j=1:length(place_fields.track)
        y_vector=[];
        for ii=1:length(place_fields.track(j).sorted_good_cells)
            %plot sorted
            matrix=[];
            normalized_matrix=[];
            matrix(ii,:)=place_fields.track(kk).rate{place_fields.track(j).sorted_good_cells(ii)};
            normalized_matrix(ii,:)=(matrix(ii,:)-min(matrix(ii,:)))/(max(matrix(ii,:))-min(matrix(ii,:)));
            subplot(length(place_fields.track),length(place_fields.track),c)
            plfield_row= normalized_matrix(ii,:)+(1.5*ii-1);
            plot(1:length(plfield_row),plfield_row,'k'); hold on;
            xx = [1:length(plfield_row), fliplr(1:length(plfield_row))];
            inBetween = [(1.5*ii-1)*ones(size(plfield_row)), fliplr(plfield_row)];
            fill(xx, inBetween,[139,0,0]/255);
            y_vector= [y_vector, 1.5*ii-1];
        end
        xlim([0 size(normalized_matrix,2)+2]);
        ylim([0 max(y_vector)+1.2]);
        yt=place_fields.track(j).sorted_good_cells;
        set(gca,'ytick',y_vector);
        set(gca,'yticklabel',yt);
        ylabel('Unit ID');
        xlabel('sorted linearized position (bins)');
        c=c+1;
        title([{['place cells on track ' num2str(j)]} ; {['sorted by track ' num2str(kk)]}]);
    end
end

% Simulate replay with linear trajectory 
ground_truth_replay = [];
sorted_spikes.replay =[];
base_duration = 0.100;
spike_times = [];
spike_id = [];
duration = [];
t_step = 0.0001; % 1ms
t = 0:t_step:1; % First 1 second is blank
spontaneous_rate=0.1;

event = 1;
for track = 1:2
    for i = 1:no_of_events
        track_id = track; % Track 1 or Track 2
        ground_truth_replay(event).track_id = track_id;
        direciton = randperm(2); % From 0 to 1 or 1 to 0
        ground_truth_replay(event).direciton = direciton(1);

        s = RandStream('mrg32k3a','Seed',track*10000+i); % Set random seed for resampling
        duration = base_duration + 2*rand(s,1)*base_duration; % Event duration ranging from 100ms to 300ms
        t_vec = t(end)+t_step:t_step:t(end)+duration;

        replay_time(1,event) = t_vec(1);
        replay_time(2,event) = t_vec(end);
        ground_truth_replay(event).duration = t_vec(end) - t_vec(1);
        ground_truth_replay(event).t = t_vec;
        t = [t t_vec(2:end)];
        
        ground_truth_replay(event).x = [];
        if ground_truth_replay(event).direciton == 1
            ground_truth_replay(event).x = linspace(0,1,length(t_vec));
        else
            ground_truth_replay(event).x = linspace(1,0,length(t_vec));
        end

        this_event_spike_times = [];
        this_event_spike_id = [];
        % Generate neuronal firing based on the 'replay trajectory'
        for j=1:length(place_fields.track(1).sorted_good_cells)
            if ground_truth_replay(event).track_id == 1
                cell = place_fields.track(1).sorted_good_cells(j);
                rate=interp1(place_fields.track(1).x,place_fields.track(1).rate{cell},ground_truth_replay(event).x);
            else
                cell = place_fields.track(2).sorted_good_cells(j);
                rate=interp1(place_fields.track(2).x,place_fields.track(2).rate{cell},ground_truth_replay(event).x);
            end

            rate=rate+spontaneous_rate;
            spike_probability=(rate*t_step).*exp(-rate*t_step);

            %     replay_spikes_max=place_field(1).track.peak_rate(j)/15;

            index=find(rand(size(t_vec))<=spike_probability);
            this_event_spike_times = [this_event_spike_times t_vec(index)];
            this_event_spike_id = [this_event_spike_id cell*ones(size(index))];

            spike_times=[spike_times t_vec(index)];
            spike_id=[spike_id cell*ones(size(index))];
        end
        
        [this_event_spike_times index] = sort(this_event_spike_times);
        ground_truth_replay(event).spikes(:,1) = this_event_spike_id(index)';
        ground_truth_replay(event).spikes(:,2) = this_event_spike_times';

        t = [t t(end):t_step:t(end)+rand(1)];
        event = event + 1;
    end
end

% Sort spike time and spike id
[spike_times index] = sort(spike_times);
replay_spikes(:,1) = spike_id(index)';
replay_spikes(:,2) = spike_times';
end
