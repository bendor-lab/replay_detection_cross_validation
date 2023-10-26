function [estimated_position,log_odds]=calculate_log_odds(spikes,place_fields,bin_size, track, opposite_track, shuffles) 

% Input
% - spikes: where spikes(:,1)
% - bins: are time bin for spike count and decoding
% - 
% By Daniel Bendor and Masahiro Takigawa 2023

common_good_cells=intersect(place_fields.track(track).sorted_good_cells,place_fields.track(opposite_track).sorted_good_cells);  %this is all common good cells
% if length(opposite_track)==1
%     common_good_cells=union(place_fields.track(track).good_cells,place_fields.track(opposite_track).good_cells);  %this is all good cells, a variant from Masa's code- less orthogonal to sequence but more sensitive
% else
%     common_good_cells=place_fields.good_place_cells;  %good place cells across all tracks
% end

% Calculate event duration
bin_edges = min(spikes(:,2)):bin_size:max(spikes(:,2));

for i=1:length(common_good_cells)
    spike_times=spikes(find(spikes(:,1)==common_good_cells(i)),2);
    binned_spikes(i,:)=histcounts(spike_times,bin_edges);
end

track_ids=[track opposite_track];
for k=1:shuffles
    for i=1:length(common_good_cells)
        s = RandStream('mrg32k3a','Seed',1000*k+i); % Set random seed for resampling

        if shuffles==1
            shuffled_track_ids=track_ids;  %don't shuffle ids if you are not calculating log odds (1 shuffle)
        else
            shuffled_track_ids=track_ids(randperm(s,length(track_ids)));
        end


        for track=1:length(track_ids)
            track_cells(track,i,:)=place_fields.track(shuffled_track_ids(track)).rate{common_good_cells(i)};  %first column is track you are analyzing (unless performing a shuffle)
        end
    end
    [estimated_position,log_odds(k)]=normalized_decoding(binned_spikes,track_cells);
   
end
 
end
    
function [estimated_position,log_odds]=normalized_decoding(binned_spikes,track_cells)
for j=1:size(binned_spikes,2)
    summed_decoded=[];
    for i=1:size(track_cells,1)
        place_fields=reshape(track_cells(i,:,:),[size(track_cells,2) size(track_cells,3)]);
        decoded(i,:)=reconstruct(binned_spikes(:,j),place_fields,0.02);
        estimated_position.track(i).decoded_positions(j,:)=decoded(i,:);
    end
    decoded=decoded/sum(sum(decoded));


    if size(track_cells,1)>2
        opposite_track(j,:)=sum(decoded(2:end,:),1);
    else
        opposite_track(j,:)=decoded(2,:);
    end

end
log_odds=log(sum(sum(estimated_position.track(1).decoded_positions))./sum(sum(opposite_track)));
end


function estimated_position = reconstruct(n,all_place_fields,bin_width)
% Creates matrix where rows are cells and columns are position bins
bin_length = size(all_place_fields,2); %columns
number_of_cells = size(all_place_fields,1); %rows
parameters.bayesian_threshold=10.^(log2(number_of_cells)-log2(400)); % small value multiplied to all values to get rid of zeros
all_place_fields(find(all_place_fields<parameters.bayesian_threshold)) = parameters.bayesian_threshold;
sum_of_place_fields = sum(all_place_fields,1);  % adds up spikes per bin (used later for exponential)
for j = 1: size(n,2)
    n_spikes = n(:,j)*ones(1,bin_length); %number of spikes in time bin
    pre_product = all_place_fields.^n_spikes; % pl field values raised to num of spikes
    pre_product(find(pre_product<parameters.bayesian_threshold)) = parameters.bayesian_threshold;
    product_of_place_fields = prod(pre_product,1); %product of pl fields
    estimated_position(:,j) = product_of_place_fields.*(exp(-bin_width*sum_of_place_fields)); % bayesian formula
    %NOTE- columns do not sum to 1.  this is done at a later stage to allow normalization within a track or across tracks
end

end