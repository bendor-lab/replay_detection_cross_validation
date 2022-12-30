function plot_confusion_matrices
% currently plots the confusion matrices for the current folder
% needs to have run the bayesian_decoding_error code first, with leave one
% out and cross tracks
est_wt= load('estimated_position_leave_one_out.mat');
est_wt= est_wt.estimated_position_test;
est_ct= load('estimated_position_cross_tracks.mat');
est_ct= est_ct.estimated_position_test;
parameters= list_of_parameters;
bin_size= parameters.x_bins_width_bayesian;


%% VERSION 1 with cross track errors
c=1;
     figure('Color','w');
     for this_track=1:length(est_wt)
         for that_track=1:length(est_wt)
            subplot(length(est_wt),length(est_wt),c);
            if this_track==that_track
                true_pos= [est_wt(this_track).run_epochs.discrete_pos_interp];
                decoded_pos= [est_wt(this_track).run_epochs.peak_position];
                n= hist3([true_pos' decoded_pos'],'Ctrs',{est_wt(this_track).position_bin_centres est_wt(this_track).position_bin_centres}); 
            else
                true_pos= [est_ct(this_track).run_epochs(that_track).true_pos_interp];
                decoded_pos= [est_ct(this_track).run_epochs(that_track).peak_position];
                n= hist3([true_pos' decoded_pos'],'Ctrs',{est_wt(this_track).position_bin_centres est_wt(this_track).position_bin_centres}); 
            end
            sum_prob= sum(n,2);
            imagesc(n./sum_prob);
            set(gca,'YDir','Normal');
            colormap gray
            map= colormap;
            colormap(flipud(map));
            title({['Track ' num2str(this_track)]; ['decoded by Track ' num2str(that_track)]});
            colorbar
            caxis([0 1]);
            xt= xticks;
            xticklabels(bin_size*xt);
            yt= yticks;
            yticklabels(bin_size*yt);
            xlabel('true position (cm)');
            ylabel('estimated position (cm)');
            c=c+1;
         end
     end
     
     
     %% VERSION 2 with normalised probabilities
     clearvars -except parameters bins_size
     load('estimated_position_leave_one_out.mat');

     figure('Color','w');
     true_pos1= [estimated_position_test(1).run_epochs(:).discrete_pos_interp];
     true_pos2= [estimated_position_test(2).run_epochs(:).discrete_pos_interp];
     true_pos= [true_pos1(1,:) true_pos2(2,:)+max(estimated_position_test(1).position_bin_centres)];
     % not ideal
     decoded_pos= [];
     for this_track=1:length(estimated_position_test)
        [~,max_indices]= max([estimated_position_test(this_track).run_epochs.max_prob]);
        peak_position_all_tmp= [estimated_position_test(this_track).run_epochs.peak_position];
        peak_position_all_tmp(2,:)= peak_position_all_tmp(2,:)+max(estimated_position_test(1).position_bin_centres);
        decoded_pos= [decoded_pos arrayfun(@(x) peak_position_all_tmp(max_indices(x),x),1:length(max_indices))];
     end
     n= hist3([true_pos' decoded_pos'],'Ctrs',{[estimated_position_test(1).position_bin_centres max(estimated_position_test(1).position_bin_centres)+estimated_position_test(1).position_bin_centres],...
                                                                            [estimated_position_test(1).position_bin_centres max(estimated_position_test(1).position_bin_centres)+estimated_position_test(1).position_bin_centres]}); 
     sum_prob= sum(n,2);
%     imagesc(log(n./sum_prob));
    imagesc(n./sum_prob);
    set(gca,'YDir','Normal');
    colormap gray
    map= colormap;
    colormap(flipud(map));
    hold on;
    plot([20.5 20.5],[0 40.5],'k','LineWidth',1.5);
    plot([0 40.5],[20.5 20.5],'k','LineWidth',1.5);
    caxis([0 1]);
    h = colorbar;
    set(get(h,'label'),'string','probability','Rotation',270);
    xlabel('true position (cm)');
    ylabel('estimated position (cm)');
    run_format_settings(gcf)
end