% SPIKE COUNT FOR BAYESIAN DECODING
% INPUTS:
% place fields: matrix or []. Place fields that want to be used as a template in the decoding. If empty, it will load place fields for the whole session (extracted_place_fields_BAYESIAN.mat)
% start_time & end_time: single value, vector or []. Might want to enter start and end timestamps when decoding only specific periods of times (e.g. laps, replay events).
%  If empty will use start and end timestamps extracted from position data (start and end times of the session).
% OUTPUTS:
% bayesian_spike_count structure. Contains time vector with time bin edges, time bin centres and matrix of spike count per place field
% Loads: 'extracted_clusters', 'extracted_place_fields_BAYESIAN','extracted_position'

function bayesian_spike_count = spike_count(place_fields_BAYESIAN,start_time,end_time,save_option)

parameters=list_of_parameters;
load extracted_clusters
load extracted_position
% If not planning to use only place fields from a specific period of time, load place fields from the whole session
if isempty(place_fields_BAYESIAN)
    load extracted_place_fields_BAYESIAN
end

%%%%%% CREATE TIME VECTOR %%%%%%
if isempty(start_time) && isempty(end_time) % if no start and end time have been input
    start_time = min(position.t);
    end_time   = max(position.t);
    t.replay_edges = start_time:parameters.replay_bin_width:end_time;     %10-20ms time bins (the smallest, the noisiest)
    t.run_edges    = start_time:parameters.run_bin_width:end_time;        %200-250ms
elseif size(start_time,2)>1 && size(end_time,2)>1 %if there's a vector for start and end times
    replay_time_edges=[];
    for i = 1: length(start_time)
        event_duration = end_time(i) - start_time(i);
        num_bins = ceil(event_duration / parameters.replay_bin_width);
        timebins_edges = linspace(start_time(i), start_time(i) + num_bins *  parameters.replay_bin_width, num_bins+1);
        replay_time_edges{i} = timebins_edges;
    end
    t.run_edges = []; % don't want to run this if just checking replay events
else  % if there's only one start and end time input
    if end_time-start_time < 1  % if decoding a replay event (short duration)
        event_duration = end_time - start_time;
        num_bins = ceil(event_duration / parameters.replay_bin_width);
        t.replay_edges = linspace(start_time, start_time + num_bins *  parameters.replay_bin_width, num_bins+1);
        t.run_edges= [];
    else
        event_duration = end_time - start_time;
        num_bins_replay = ceil(event_duration / parameters.replay_bin_width);
        num_bins_run = ceil(event_duration / parameters.run_bin_width);
        t.replay_edges = linspace(start_time, start_time + num_bins_replay *  parameters.replay_bin_width, num_bins_replay+1);
        t.run_edges    = linspace(start_time, start_time + num_bins_run *  parameters.run_bin_width, num_bins_run+1);
    end
end

%select place cells to use for replay analysis

if isfield(place_fields_BAYESIAN,'good_place_cells')
    place_field_index = place_fields_BAYESIAN.good_place_cells;  %use all place cells
elseif ~isfield(place_fields_BAYESIAN,'good_place_cells') && isfield(place_fields_BAYESIAN.track,'good_cells')
    place_field_index = place_fields_BAYESIAN.track.good_cells;  % when using output from get_place_fields_laps
else
    disp('ERROR- field good_place_cells missing');
    place_field_index = place_fields_BAYESIAN.pyramidal_cells;
end

%%%%% Create spikes count matrix %%%%%
bayesian_spike_count.n.replay= [];
if exist('replay_time_edges','var')  % When running replay events separately
    bayesian_spike_count.replay_events_indices =[];   replay_edges_concat=[];
    for i = 1 : length(replay_time_edges)
        t.replay_edges = cell2mat(replay_time_edges(i)); % takes each replay time vector separately
        % Takes time vectors and centres each bin
        bayesian_spike_count.replay_events(i).replay_time_edges = t.replay_edges;
        replay_edges_concat=[replay_edges_concat t.replay_edges];
        bayesian_spike_count.replay_events(i).replay_time_centered = t.replay_edges(1:end-1)+parameters.replay_bin_width/2; %centres of bins
        %concatenate indicies for each replay event, put NaN in between events so that spikes between edges of replay events are ignored
        bayesian_spike_count.replay_events_indices = [bayesian_spike_count.replay_events_indices, ones(1,length(bayesian_spike_count.replay_events(i).replay_time_centered))*i,NaN];
    end
    % Spike histogram per time bin for each place field,
    % performed on bins across all replay events (time bin between replay events gets ignored later because replay_event_indice is set to NaN
    for k=1:length(place_field_index)
        bayesian_spike_count.n.replay(k,:) = histcounts(clusters.spike_times(find(clusters.spike_id==place_field_index(k))),replay_edges_concat);
    end
else   % When running whole sesion together- Takes time vectors and centres each bin
    bayesian_spike_count.replay_time_edges = t.replay_edges;
    bayesian_spike_count.replay_time_centered = t.replay_edges(1:end-1)+parameters.replay_bin_width/2; %centres of bins
    % Spike histogram per time bin for each place field (rows = cells /columns = time bin)
    for k=1:length(place_field_index)
        bayesian_spike_count.n.replay(k,:) = histcounts(clusters.spike_times(find(clusters.spike_id==place_field_index(k))),t.replay_edges);
    end
    %compute also for run if needed
    if ~isempty(t.run_edges)
        bayesian_spike_count.run_time_edges = t.run_edges;
        bayesian_spike_count.run_time_centered = t.run_edges(1:end-1)+parameters.run_bin_width/2;%centres of bins
        bayesian_spike_count.n.run = [];
        for k=1:length(place_field_index)
            bayesian_spike_count.n.run(k,:) = histcounts(clusters.spike_times(find(clusters.spike_id==place_field_index(k))),t.run_edges);
        end
    end 
end

if strcmp(save_option, 'Y')
    save bayesian_spike_count bayesian_spike_count
end
end