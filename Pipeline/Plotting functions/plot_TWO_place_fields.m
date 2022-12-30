function plot_TWO_place_fields(place_fields1,place_fields2,track_id)
place_fields{1}=place_fields1.track(track_id);
place_fields{2}=place_fields2.track(track_id);
%%% PLOT WITH IMAGESC
c=1;
figure
y_vector=[];
for i=1:2
    for j=1:2
        matrix=[];
        normalized_matrix=[];
        subplot(2,2,c)
        for k=1:length(place_fields{i}.sorted_good_cells)
            matrix(k,:)=place_fields{i}.smooth{place_fields{j}.sorted_good_cells(k)};
            normalized_matrix(k,:)=(matrix(k,:)-min(matrix(k,:)))/(max(matrix(k,:))-min(matrix(k,:)));
        end
        y_vector= [y_vector, 1.5*i-1];
        imagesc(normalized_matrix);
        colormap(jet)
        hold on
        ylabel('cell id');
        yt=place_fields{j}.sorted_good_cells;
        set(gca,'yticklabel',yt);
        xlabel('linearized position')
        title([{['place cells on track ' num2str(j)]} ; {['sorted by track ' num2str(i)]}]);
        axis xy
        c=c+1;
    end
end

figure
plot(place_fields{1}.raw_peak(place_fields{1}.good_cells),place_fields{2}.raw_peak(place_fields{2}.good_cells),'o')
[r,p]=corr(place_fields{1}.raw_peak(place_fields{1}.good_cells)',place_fields{2}.raw_peak(place_fields{2}.good_cells)','type','Spearman')
end