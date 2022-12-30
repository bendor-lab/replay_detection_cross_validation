function plot_half_session_place_fields()
% helps create plots, with the ratemaps divided into first and second half
% (use laps)
% to show location is stable

parameters= list_of_parameters;
load('lap_times.mat');
load('extracted_place_fields_BAYESIAN.mat');

% find cells common to both tracks
[~,idx_T1,idx_T2]= intersect(place_fields_BAYESIAN.track(1).good_cells,place_fields_BAYESIAN.track(2).good_cells);
cdata= [{repmat([0.5 0.5 0.5],length(place_fields_BAYESIAN.track(1).good_cells),1)} {repmat([0.5 0.5 0.5],length(place_fields_BAYESIAN.track(2).good_cells),1)}];
cdata{1}(idx_T1,3)= 1; cdata{1}(idx_T1,1:2)= 0;
cdata{2}(idx_T2,3)= 1; cdata{2}(idx_T2,1:2)= 0;

for this_track=1:length(lap_times)
    middle_lap= round(lap_times(this_track).total_number_of_laps/2); % roughly 50% of laps
    % find starts & stops 
    h1{this_track,:}= [lap_times(this_track).start(1) lap_times(this_track).end(middle_lap-1) ];
    h2{this_track,:}= [lap_times(this_track).start(middle_lap) lap_times(this_track).end(end)];
end

    % calculate place fields for each half
    pl_fld_h1= calculate_place_fields_epochs(parameters.x_bins_width_bayesian,h1);
    pl_fld_h2= calculate_place_fields_epochs(parameters.x_bins_width_bayesian,h2);
%     
%     % plot
%     c=1;
%     fd= pwd;
%     fd= strsplit(fd,'\');
%     fd= [fd{end-1} '\' fd{end}];
%     figure('Name',fd,'Color',[1 1 1])
%     y_vector=[];
%     for i=1:length(place_fields_BAYESIAN.track)
%         for j=1:length(place_fields_BAYESIAN.track)
%             matrix=[];
%             normalized_matrix=[];
%             subplot(length(place_fields_BAYESIAN.track),length(place_fields_BAYESIAN.track),c)
%             for k=1:length(place_fields_BAYESIAN.track(j).sorted_good_cells)
%                 % collate halves
%                 matrix(k,:)= [pl_fld_h1.track(i).raw{place_fields_BAYESIAN.track(j).sorted_good_cells(k)} pl_fld_h2.track(i).raw{place_fields_BAYESIAN.track(j).sorted_good_cells(k)}];
%                 normalized_matrix(k,:)=(matrix(k,:)-min(matrix(k,:)))/(max(matrix(k,:))-min(matrix(k,:)));
%             end
%             y_vector= [y_vector, 1.5*i-1];
%             imagesc(normalized_matrix);
% %             colormap(jet)
%             colormap(pink)
%             hold on;
%             plot([size(matrix,2)/2+0.5 size(matrix,2)/2+0.5],ylim,'w--','LineWidth',2);
%             ylabel('cell id');
%             yticks(1:5:length(place_fields_BAYESIAN.track(j).sorted_good_cells));
%             yticklabels(1:5:length(place_fields_BAYESIAN.track(j).sorted_good_cells));
%             if c>6
%                 xlabel('linearized position (cm)')
%             end
%             xticks([1:4:size(matrix,2)]);
%             xticklabels([pl_fld_h1.track(this_track).x_bin_edges(1:4:end)])
%             title([{['place cells on track ' num2str(j)]} ; {['sorted by track ' num2str(i)]} ; {'   first half         second half   '}]);
%             axis xy
%             c=c+1;
%         end
%     end
%     
    c=1;
    figure('Color','w');
     for i=1:length(place_fields_BAYESIAN.track)
        for j=1:length(place_fields_BAYESIAN.track)
            matrix=[];
            normalized_matrix=[];
            subplot(length(place_fields_BAYESIAN.track),length(place_fields_BAYESIAN.track),c)
            for k=1:length(place_fields_BAYESIAN.track(j).sorted_good_cells)
                % collate halves
                matrix(k,:)= [pl_fld_h1.track(i).raw{place_fields_BAYESIAN.track(j).sorted_good_cells(k)} pl_fld_h2.track(i).raw{place_fields_BAYESIAN.track(j).sorted_good_cells(k)}];
                normalized_matrix(k,:)=(matrix(k,:)-min(matrix(k,:)))/(max(matrix(k,:))-min(matrix(k,:)));
            end
            hold on;
            arrayfun(@(x) plot(x+ normalized_matrix(x,:),'Color',cdata{j}(x,:)),1:size(normalized_matrix,1));
            arrayfun(@(x) fill([1:size(normalized_matrix,2) fliplr(1:size(normalized_matrix,2))],[x*ones(1,size(normalized_matrix,2)) fliplr(x+normalized_matrix(x,:))],cdata{j}(x,:)),1:size(normalized_matrix,1));
            plot([size(matrix,2)/2+0.5 size(matrix,2)/2+0.5],ylim,'w--','LineWidth',2);
             ylim([1 size(normalized_matrix,1)+1])
            if c==1
                xlabel('linearized position (cm)')
            end
            if ismember(c,[1 4])
                ylabel(['number of cells: ' num2str(length(place_fields_BAYESIAN.track(j).sorted_good_cells))]);
                yticks([]);
                    set(gca,'TickDir','out','XColor','none');
                  if c==1
                    xticks([1 size(matrix,2)/2]);
                    xticklabels([pl_fld_h1.track(this_track).x_bin_edges([1 end])])
                    set(gca,'TickDir','out','XColor','k');
                  end
            else 
                set(gca,'Xcolor','none','YColor','none');
            end
            if c==1
                title([{['units on Track ' num2str(j)]};{['sorted by Track ' num2str(i)]} ; {'   first 50% laps                             remaining 50% laps   '}]);
            else
                title([{['units on Track ' num2str(j)]};{['sorted by Track ' num2str(i)]}]);
            end
            axis xy
            c=c+1;
        end
     end
    



end