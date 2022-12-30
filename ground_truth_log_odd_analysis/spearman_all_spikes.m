function [corr_value, p_value]=spearman_all_spikes(spike_id, spike_times, sorted_place_fields)

 cell_IDs = zeros(1,length(spike_id));
 
 for cell = 1 : length(sorted_place_fields)
     index = find(spike_id == sorted_place_fields(cell));
     if ~isempty(index)
         cell_IDs(index) = cell;
     end
 end

 if isempty(cell_IDs)
     corr_value  = NaN;
     p_value = NaN;
 else
     
     [corr_value, p_value] = corr(cell_IDs',spike_times,'type','Spearman');
     corr_value = abs(corr_value);  %just care about magnitude of correlation
 end

end