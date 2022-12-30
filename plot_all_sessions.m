function plot_all_sessions(varargin)

cd 'X:\BendorLab\Drobo\Neural and Behavioural Data\Rate remapping\Data'
load folders_to_process_remapping

parameters = list_of_parameters;
f1 = figure;
    data_across_tracks(1).place_field_diff=[];
    data_across_tracks(1).mean_spike_diff=[];
    data_across_tracks(1).mean_spike_diff_nonZero=[];
    data_across_tracks(1).track_id=[];
    data_across_tracks(2).place_field_diff=[];
    data_across_tracks(2).mean_spike_diff=[];
    data_across_tracks(2).mean_spike_diff_nonZero=[];
    data_across_tracks(2).track_id=[];
    
for i = 1 : length(folders)
    
    cd(['X:\BendorLab\Drobo\Neural and Behavioural Data\Rate remapping\Data\' folders{i}])
        
    % Load track data
    if ~isempty(varargin)
        switch varargin{1}
            case 'wcorr'
                load  rate_remapping_analysis_TRACK_PAIRS_wcorr
            case 'spearman'
                load rate_remapping_analysis_TRACK_PAIRS_spearman
        end
    end
    

    %%%%% PLOT T1
    figure(f1)

    for epoch=1:size(remapping,1)
        
        for track_pair = 1:size(remapping,2)
            [~,common_PRE_Post_indx,~] = intersect(remapping(epoch,track_pair).ID_active_cells_during_replay,remapping(epoch,track_pair).PRE_to_POST_active_cells);
            subplot(size(remapping,1),2,epoch)
            hold on
            plot(remapping(epoch,track_pair).place_field_diff, remapping(epoch,track_pair).replay_spike_diff_nonZero,'o','MarkerEdgeColor',[0.2 0.4 0.8]);
            hold on
            plot(remapping(epoch,track_pair).place_field_diff(common_PRE_Post_indx), remapping(epoch,track_pair).replay_spike_diff_nonZero(common_PRE_Post_indx),'o','MarkerEdgeColor',[0.2 0.4 0.8],'MarkerFaceColor',[0.2 0.4 0.8]);
            subplot(size(remapping,1),2,epoch+2)
            plot(remapping(epoch,track_pair).place_field_diff, remapping(epoch,track_pair).replay_spike_diff,'o','MarkerEdgeColor',[0.2 0.4 0.8]);
            hold on
            plot(remapping(epoch,track_pair).place_field_diff(common_PRE_Post_indx), remapping(epoch,track_pair).replay_spike_diff(common_PRE_Post_indx),'o','MarkerEdgeColor',[0.2 0.4 0.8],'MarkerFaceColor',[0.2 0.4 0.8]);
            
            % index=find(abs(remapping(epoch,track_pair).place_field_centre_diff)<=10);  %place field centres less than 20 cm apart
            % plot(remapping(epoch,track_pair).place_field_diff(index), remapping(epoch,track_pair).replay_spike_diff(index),'k.');
            index_non_NaNs=find(~isnan(remapping(epoch,track_pair).place_field_diff));
            data_across_tracks(epoch).place_field_diff=[data_across_tracks(epoch).place_field_diff; remapping(epoch,track_pair).place_field_diff(index_non_NaNs)];
            data_across_tracks(epoch).mean_spike_diff=[data_across_tracks(epoch).mean_spike_diff; remapping(epoch,track_pair).replay_spike_diff(index_non_NaNs)];
            data_across_tracks(epoch).mean_spike_diff_nonZero=[data_across_tracks(epoch).mean_spike_diff_nonZero; remapping(epoch,track_pair).replay_spike_diff_nonZero(index_non_NaNs)];
            
            %[r,p]=corr(remapping(epoch,track_pair).place_field_diff(index_non_NaNs), remapping(epoch,track_pair).replay_spike_diff(index_non_NaNs),'type','Pearson');
            %Epoch_and_RandP_values=[epoch,track_pair,r]
            %p
        end
        if i ==5
            if ~isempty(data_across_tracks(epoch).place_field_diff) & ~isempty(data_across_tracks(epoch).mean_spike_diff)
                lm = fitlm(data_across_tracks(epoch).place_field_diff, data_across_tracks(epoch).mean_spike_diff_nonZero,'linear')
                [p,~,~] = coefTest(lm);
                linear_fit_P_values_RATE_VS_PEAK_PF_nonZero(epoch) = p;
                x=[min(data_across_tracks(epoch).place_field_diff) max(data_across_tracks(epoch).place_field_diff)];
                b=lm.Coefficients.Estimate';
                subplot(size(remapping,1),2,epoch)
                plot(x,polyval(fliplr(b),x),'k--');
                
                
                
                lm = fitlm(data_across_tracks(epoch).place_field_diff, data_across_tracks(epoch).mean_spike_diff,'linear')
                [p,~,~] = coefTest(lm);
                linear_fit_P_values_RATE_VS_PEAK_PF(epoch)= p;
                x=[min(data_across_tracks(epoch).place_field_diff) max(data_across_tracks(epoch).place_field_diff)];
                b=lm.Coefficients.Estimate';
                subplot(size(remapping,1),2,epoch+2)
                plot(x,polyval(fliplr(b),x),'k--');
            end
        end
       
    end
end

    ax1 = subplot(2,2,1)
    title(['sleep- PRE nonzero , pval= ' num2str(linear_fit_P_values_RATE_VS_PEAK_PF_nonZero(1))])
    ylabel({'Mean # spikes diff during replay'});
    xlabel('diff place field peak FR')
    box off
    ax1.FontSize = 16;

    ax2 = subplot(2,2,2)
    title(['sleep- POST nonzero , pval= ' num2str(linear_fit_P_values_RATE_VS_PEAK_PF_nonZero(2))])
    ylabel({'Mean # spikes diff during replay'});
    xlabel('diff place field peak FR')
    box off
    ax2.FontSize = 16;

    ax3 = subplot(2,2,3)
    title(['sleep- PRE, pval= ' num2str(linear_fit_P_values_RATE_VS_PEAK_PF(1))])
     ylabel({'Mean # spikes diff during replay'});
    xlabel('diff place field peak FR')
    box off
    ax3.FontSize = 16;

    ax4 = subplot(2,2,4)
    title(['sleep- POST, pval= ' num2str(linear_fit_P_values_RATE_VS_PEAK_PF(2))])
     ylabel({'Mean # spikes diff during replay'});
    xlabel('diff place field peak FR')
    box off
    ax4.FontSize = 16;
 end

