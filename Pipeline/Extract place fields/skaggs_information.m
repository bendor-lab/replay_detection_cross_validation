function skaggs_info= skaggs_information(varargin)
% MT - adapted from Skaggs, McNaughton, Gothard
% calculates skaggs information according to following formula:
% info= sum_over_all_x(mean_firing_rate_at_x*log2(mean_firing_rate_at_x/overall_mean_firing_rate_of_cell)*prob_density_being_at_x)

if isempty(varargin)
    load extracted_place_fields.mat
    for i=1:length(place_fields.track)
        for j=1:length(place_fields.track(i).sorted_good_cells)

            %go through each good sorted cell
            cell= place_fields.track(i).sorted_good_cells(j);
            %need rate map
            ratemap= cell2mat(place_fields.track(i).raw(cell));
            %need dwell map
            dwellmap= place_fields.track(i).dwell_map;

            %remove locations where the firing rate is 0
            dwellmap(ratemap==0)=[];
            ratemap(ratemap==0)=[];

            % probability of being at location x
            prob_x= dwellmap/nansum(dwellmap);

            % overall firing rate on track
            meanFR= nansum(ratemap.*prob_x);

            norm_rate= ratemap/meanFR;
            log_norm= log2(norm_rate);

            place_fields.track(i).skaggs_info_sorted_good_cells(j)= nansum(prob_x.*norm_rate.*log_norm);
        end
    end

    save('extracted_place_fields.mat','place_fields');
    
elseif length(varargin)==1
    
    fields= varargin(1);
    fields= fields{1,1}; %too lazy to fix it better
    for j=1:length(fields.raw)
            
            %need rate map
            ratemap= cell2mat(fields.raw(j));
            %need dwell map
            dwellmap= fields.dwell_map;

            %remove locations where the firing rate is 0
            dwellmap(ratemap==0)=[];
            ratemap(ratemap==0)=[];

            % probability of being at location x
            prob_x= dwellmap/nansum(dwellmap);

            % overall firing rate on track
            meanFR= nansum(ratemap.*prob_x);
            
%             if meanFR ~= fields.mean_rate_track(j)
%                 keyboard;
%             end

            norm_rate= ratemap/meanFR;
            log_norm= log2(norm_rate);

            skaggs_info(j)= nansum(prob_x.*norm_rate.*log_norm);
    end
end

end