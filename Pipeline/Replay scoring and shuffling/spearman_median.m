function [corr_value, p_value]=spearman_median(spike_id, spike_times, sorted_place_fields)
 c = 1;
 cell_IDs = []; median_spikeTime = [];
 
 for cell = 1 : length(sorted_place_fields)
     index = find(spike_id == sorted_place_fields(cell));
     if ~isempty(index)
         cell_IDs(c) = cell;
         median_spikeTime(c) = median(spike_times(index));
         c = c+1;
     end
 end

 if isempty(cell_IDs) || isempty(median_spikeTime)
     corr_value  = NaN;
     p_value = NaN;
 else
     
 [corr_value, p_value] = corr(cell_IDs',median_spikeTime','type','Spearman');
 corr_value = abs(corr_value);  %just care about magnitude of correlation
 end

    