function place_fields = calculate_place_fields(x_bins_width)
% INPUTS:
%   x_bin_width: enter width value (2 for fine resolution, 10 for bayesian
%   decoding).
% loads list_of_parameters.m, extracted_clusters.mat,extracted_position.mat, extracted_waveform.mat
% uses function skaggs_information.m
% Main difference is to use bin centres rather than edges

    parameters = list_of_parameters;
    position = [];
    clusters = [];
    load('extracted_clusters.mat');
    load('extracted_position.mat');
    if exist('extracted_waveforms.mat','file')
        load('extracted_waveforms.mat');
    else
        disp('no extracted_waveforms.mat file');
        allclusters_waveform=[];
    end

    %find track positions
    for track_id=1:length(position.linear)
        track_indices = ~isnan(position.linear(track_id).linear);
        place_fields.track(track_id).time_window=[min(position.t(track_indices)) max((position.t(track_indices)))];
    end

    % Run threshold on pyramidal cells: half-width amplitude
    if ~isempty(allclusters_waveform)       
        PC_indices = [allclusters_waveform.half_width] > parameters.half_width_threshold; % cells that pass treshold of pyramidal cell half width
        pyramidal_cells = [allclusters_waveform(PC_indices).converted_ID];
    end
    
    %find mean rate spikes on all tracks / time on all tracks
    %(for identifying putative interneurons later in the code)
    for track_id = 1:length(position.linear)
        total_time_in_track(track_id) = place_fields.track(track_id).time_window(2)-place_fields.track(track_id).time_window(1);
        for j = 1 : max(clusters.id_conversion(:,1))
            all_spikes(track_id,j) = length(find(clusters.spike_id==j & ...
                clusters.spike_times>place_fields.track(track_id).time_window(1) & ...
                clusters.spike_times<place_fields.track(track_id).time_window(2)));
        end
    end
    place_fields.mean_rate=sum(all_spikes,1)/sum(total_time_in_track);
    
    
    %% Place field calculation
    
    time_bin_width=position.t(2)-position.t(1);

    for track_id = 1:length(position.linear)
        
        position_index = isnan(position.linear(track_id).linear);
        position_speed = abs(position.v_cm);
        position_speed(position_index) = NaN;  %make sure speed is NaN if position is NaN

        position_during_spike = interp1(position.t,position.linear(track_id).linear,clusters.spike_times,'nearest'); %interpolates position into spike time
        speed_during_spike = interp1(position.t,position_speed,clusters.spike_times,'nearest');

        x_bin_edges = 0:x_bins_width:100*position.linear(track_id).length; % forces x_bins to be from 0 to 200cm
        x_bin_centres = [(x_bin_edges(2)-x_bins_width/2):x_bins_width:(x_bin_edges(end-1)+x_bins_width/2)];
        x_bins = 0:x_bins_width:(100*position.linear(track_id).length); %bin position

        % Time spent at each x_bin (speed filtered)
        x_hist = time_bin_width.*histcounts(position.linear(track_id).linear(find(position.t>place_fields.track(track_id).time_window(1) &...
                 position.t<place_fields.track(track_id).time_window(2) & position_speed>parameters.speed_threshold_laps...
                 & position_speed<parameters.speed_threshold_max)),x_bin_edges); % Changed bin_centre to bin_edges

        place_fields.track(track_id).x_bin_centres = x_bin_centres;
        place_fields.track(track_id).x_bin_edges = x_bin_edges;
        place_fields.track(track_id).x_bins = x_bins;
        place_fields.track(track_id).x_bins_width = x_bins_width;
        place_fields.track(track_id).dwell_map = x_hist;

        for j = 1 : max(clusters.id_conversion(:,1))
            
            % Number of spikes per bin within time window (speed filtered)
            place_fields.track(track_id).spike_hist{j} = histcounts(position_during_spike(find(clusters.spike_id==j & ...
                                                        clusters.spike_times>place_fields.track(track_id).time_window(1) & ...
                                                        clusters.spike_times<place_fields.track(track_id).time_window(2) & ...
                                                        speed_during_spike>parameters.speed_threshold_laps &...
                                                        speed_during_spike<parameters.speed_threshold_max)),x_bin_edges); % Changed bin_centre to bin_edges

            place_fields.track(track_id).raw{j} = place_fields.track(track_id).spike_hist{j}./x_hist; % place field calculation
            place_fields.track(track_id).raw{j}(find(isnan(place_fields.track(track_id).raw{j})==1))=0;
          
            % zero bins with 0 dwell time, but make sure no spikes occurred
            non_visited_bins = find(x_hist==0);
            if sum(place_fields.track(track_id).spike_hist{j}(non_visited_bins))>0
                disp('ERROR: x_hist is zero, but spike histogram is not');
            else
                place_fields.track(track_id).raw{j}(non_visited_bins)= 0;
            end
            place_fields.track(track_id).non_visited_bins = non_visited_bins; %NaNs that have been replaced by O
            
            % Create smoothing filter (gamma)
            if x_bins_width== parameters.x_bins_width_bayesian
                w= [1 1];  %moving average filter of 2 sample, will be become a filter of [0.25 0.5 0.25] with filtfilt
            else
                w= gausswin(parameters.place_field_smoothing);
            end
            w = w./sum(w); %make sure smoothing filter sums to 1
            
            % Get place field information
            place_fields.track(track_id).smooth{j}         = filtfilt(w,1,place_fields.track(track_id).raw{j}); %smooth pl field
            place_fields.track(track_id).centre_of_mass(j) = sum(place_fields.track(track_id).smooth{j}.*x_bin_centres/sum(place_fields.track(track_id).smooth{j}));  %averaged center
            [place_fields.track(track_id).peak(j) , index] = max(place_fields.track(track_id).smooth{j}); %peak of smoothed place field and index of peak (center)
            if place_fields.track(track_id).peak(j) ~=0
                if length(index)>1 % very rare exception where you have multiple peaks of same height....
                    index= index(1);
                end
                place_fields.track(track_id).centre(j) = x_bin_centres(index);
                
            else
                place_fields.track(track_id).centre(j) = NaN;
            end
            place_fields.track(track_id).raw_peak(j)          = max(place_fields.track(track_id).raw{j}); % raw pl field peak
            place_fields.track(track_id).mean_rate_session(j) = length(find(clusters.spike_id==j))/(position.t(end)-position.t(1)); %mean firing rate
            place_fields.track(track_id).mean_rate_track(j)   = sum(place_fields.track(track_id).spike_hist{j})/(place_fields.track(track_id).time_window(2)-place_fields.track(track_id).time_window(1));
            if place_fields.track(track_id).peak(j) ~=0
                place_fields.track(track_id).half_max_width(j) = x_bins_width*half_max_width(place_fields.track(track_id).smooth{j}); %finds half width of smoothed place field (width from y values closest to 50% of peak)
            else
                place_fields.track(track_id).half_max_width(j) = NaN;
            end
        end

        %calculate skagges information
        place_fields.track(track_id).skaggs_info= skaggs_information(place_fields.track(track_id));
        
        % Find cells that pass the 'Place cell' thresholds -
        % both peak of smoothed place field or peak of raw place field need to be above the respective thresholds
        putative_place_cells = find((place_fields.track(track_id).peak >= parameters.min_smooth_peak...
            & place_fields.track(track_id).raw_peak >= parameters.min_raw_peak)...
            & place_fields.mean_rate <= parameters.max_mean_rate...
            & place_fields.track(track_id).skaggs_info > 0);
        
        % Set a less conservative criteria for place cells, having to pass either peak firing rate thresholds (smoothed PF and raw PF)
        putative_place_cells_LIBERAL = find(place_fields.track(track_id).peak >= parameters.min_smooth_peak... 
            | place_fields.track(track_id).raw_peak >= parameters.min_raw_peak...
            & place_fields.mean_rate <= parameters.max_mean_rate...
            & place_fields.track(track_id).skaggs_info > 0);
        
        if ~isempty(allclusters_waveform)
            place_fields.track(track_id).good_cells = intersect(putative_place_cells,pyramidal_cells); % Check that the cells that passed the threshold are pyramidal cells
            place_fields.track(track_id).good_cells_LIBERAL = intersect(putative_place_cells_LIBERAL,pyramidal_cells); % Check that the cells that passed the threshold are pyramidal cells
        else
            place_fields.track(track_id).good_cells = putative_place_cells;
            place_fields.track(track_id).good_cells_LIBERAL = putative_place_cells_LIBERAL;            
        end
       
        % Sort place fields according to the location of their peak
        [~,index] = sort(place_fields.track(track_id).centre);
        place_fields.track(track_id).sorted = index;
        [~,index1] = sort(place_fields.track(track_id).centre(place_fields.track(track_id).good_cells));
        place_fields.track(track_id).sorted_good_cells = place_fields.track(track_id).good_cells(index1);
        [~,index2] = sort(place_fields.track(track_id).centre(place_fields.track(track_id).good_cells_LIBERAL));
        place_fields.track(track_id).sorted_good_cells_LIBERAL = place_fields.track(track_id).good_cells_LIBERAL(index2);
    end
    
    %% Classify cells as good place cells, interneuron, pyramidal cells & other cells
    
     %interneurons classfication
     interneurons = find(place_fields.mean_rate > parameters.max_mean_rate);
     place_fields.interneurons=interneurons;
              
    good_place_cells=[]; track=[];
    for track_id=1:length(position.linear) %good cells classfication
        all_cells = 1:1:length(place_fields.track(track_id).raw);
        good_place_cells = [good_place_cells place_fields.track(track_id).sorted_good_cells];
        track =[track track_id*ones(size(place_fields.track(track_id).sorted_good_cells))];
    end
    place_fields.good_place_cells = unique(good_place_cells);
    place_fields.all_cells = unique(all_cells);
    
    good_place_cells_LIBERAL=[];
    for track_id=1:length(position.linear) %good cells (liberal threshold) classfication
        good_place_cells_LIBERAL = [good_place_cells_LIBERAL place_fields.track(track_id).sorted_good_cells_LIBERAL];
    end
    place_fields.good_place_cells_LIBERAL = unique(good_place_cells_LIBERAL);
    
    % cells that are unique for each track
    unique_cells=[];
    for track_id = 1:length(position.linear)
        place_fields.track(track_id).unique_cells = setdiff(good_place_cells(track==track_id),good_place_cells(track~=track_id),'stable');
        unique_cells = [unique_cells, place_fields.track(track_id).unique_cells];
    end
    place_fields.unique_cells = unique_cells;  % all cells that have good place fields only on a single track
    
   % putative pyramidal cells classification:  pyramidal cells that pass the 'Pyramidal type' threshold (but not need to be place cells)
    putative_pyramidal_cells = find(place_fields.mean_rate <= parameters.max_mean_rate);
    
    if ~isempty(allclusters_waveform)
        place_fields.pyramidal_cells = intersect(putative_pyramidal_cells,pyramidal_cells);
    else
        place_fields.pyramidal_cells = putative_pyramidal_cells;
    end
    place_fields.pyramidal_cells=unique(place_fields.pyramidal_cells);

    other_cells = setdiff(1:max(clusters.id_conversion(:,1)),good_place_cells,'stable'); %find the excluded putative pyramidal cells
    place_fields.other_cells = setdiff(other_cells,interneurons,'stable'); %remove also the interneurons
    
    %save place fields (different filenames used based on x_bins_width chosen)
    if x_bins_width== parameters.x_bins_width_bayesian
        place_fields_BAYESIAN=place_fields;
        save extracted_place_fields_BAYESIAN place_fields_BAYESIAN;
    elseif x_bins_width== parameters.x_bins_width
        save extracted_place_fields place_fields;
    else disp('error: x_bin_width does not match expected value')
    end
    
end

function half_width = half_max_width(place_field)
    %interpolate place field to get better resolution
    new_step_size = 0.1;  %decrease value to get finer resolution interpolation of place field
    place_field_resampled = interp1(1:length(place_field),place_field,1:new_step_size:length(place_field),'linear');
    [peak,index] = max(place_field_resampled); %finds smoothed place field peak firing rate (FR)
    for i = index : length(place_field_resampled)
        if place_field_resampled(i)<peak/2 %finds the point after the peak where the FR is half the peak FR
            break;
        end
    end
    for j = index : -1 : 1 %finds the point before the peak where the FR is half the peak FR
        if place_field_resampled(j)<peak/2
            break;
        end
    end
    half_width = new_step_size*(i-j); %distance between half-peaks
    %(calculated in indicies of original place field, but converted to distance in cm in function above)
end
