function plot_rate_remapping(option,varargin)
f1= figure('units','normalized','outerposition',[0 0 1 1]);
parameters=list_of_parameters;
if strcmp(option,'TRACK_PAIRS')
    
    if ~isempty(varargin)
        switch varargin{1}
            case 'wcorr'
                load('rate_remapping_analysis_TRACK_PAIRS_wcorr');
            case 'spearman'
                load('rate_remapping_analysis_TRACK_PAIRS_spearman');
            otherwise
                load('rate_remapping_analysis_TRACK_PAIRS');
        end
    else
        varargin{1}=[];
        load('rate_remapping_analysis_TRACK_PAIRS');
    end
    
    foldername = strsplit(pwd,'\');
    f1.Name= strcat(varargin{1},'_',foldername{7});
    
    for epoch=1:size(remapping,1)
        data_across_tracks(epoch).place_field_diff=[];
        data_across_tracks(epoch).mean_spike_diff=[];
        data_across_tracks(epoch).mean_spike_diff_nonZero=[];
        data_across_tracks(epoch).track_id=[];
        for track_pair = 1:size(remapping,2)
            [~,common_PRE_Post_indx,~] = intersect(remapping(epoch,track_pair).ID_active_cells_during_replay,remapping(epoch,track_pair).PRE_to_POST_active_cells);
            figure(f1)
            subplot(size(remapping,1),2,epoch)
            plot(remapping(epoch,track_pair).place_field_diff, remapping(epoch,track_pair).replay_spike_diff_nonZero,parameters.plot_color_symbol{track_pair});
            hold on
            plot(remapping(epoch,track_pair).place_field_diff(common_PRE_Post_indx), remapping(epoch,track_pair).replay_spike_diff_nonZero(common_PRE_Post_indx),'k*');
            subplot(size(remapping,1),2,epoch+2)
            plot(remapping(epoch,track_pair).place_field_diff, remapping(epoch,track_pair).replay_spike_diff,parameters.plot_color_symbol{track_pair});
            hold on
            plot(remapping(epoch,track_pair).place_field_diff(common_PRE_Post_indx), remapping(epoch,track_pair).replay_spike_diff(common_PRE_Post_indx),'k*');
            
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
    
    
    subplot(2,2,1)
    title(['sleep- PRE nonzero , pval= ' num2str(linear_fit_P_values_RATE_VS_PEAK_PF_nonZero(1))])
    ylabel([{'Average spike number difference during replay between tracks'}; {'divided by number of events cell was active in)'}])
    xlabel('place field peak response difference between tracks')
    subplot(2,2,2)
    title(['sleep- POST nonzero , pval= ' num2str(linear_fit_P_values_RATE_VS_PEAK_PF_nonZero(2))])
    ylabel([{'Average spike number difference during replay between tracks'}; {'divided by number of events cell was active in)'}])
    xlabel('place field peak response difference between tracks')
    subplot(2,2,3)
    title(['sleep- PRE, pval= ' num2str(linear_fit_P_values_RATE_VS_PEAK_PF(1))])
    ylabel('Average spike number difference during replay between tracks')
    xlabel('place field peak response difference between tracks')
    subplot(2,2,4)
    title(['sleep- POST, pval= ' num2str(linear_fit_P_values_RATE_VS_PEAK_PF(2))])
    ylabel('Average spike number difference during replay between tracks')
    xlabel('place field peak response difference between tracks')
end


if strcmp(option,'ONE_TRACK')
    
     if ~isempty(varargin)
        switch varargin{1}
            case 'wcorr'
                load('rate_remapping_analysis_ONE_TRACK_wcorr');
            case 'spearman'
                load('rate_remapping_analysis_ONE_TRACK_spearman');
            otherwise
                load('rate_remapping_analysis_ONE_TRACK');
        end
    else
        load('rate_remapping_analysis_ONE_TRACK');
     end
    
    for epoch=1:size(remapping,1)
        data_across_tracks(epoch).place_field_peak=[];
        data_across_tracks(epoch).mean_replay_rate=[];
        data_across_tracks(epoch).track_id=[];
        for track=1:size(remapping,2)
            
            
            index_active_place_fields=find(remapping(epoch,track).place_fields_BAYESIAN>=1);  %only analyze neurons with place fields on tracks
            index_NOTactive_place_fields_that_replay=find(remapping(epoch,track).place_fields_BAYESIAN==0 & remapping(epoch,track).mean_replay_rate'>0);
            if ~isempty(index_active_place_fields)
                data_across_tracks(epoch).place_field_peak=[data_across_tracks(epoch).place_field_peak remapping(epoch,track).place_fields_BAYESIAN(index_active_place_fields)];
                data_across_tracks(epoch).mean_replay_rate=[data_across_tracks(epoch).mean_replay_rate; remapping(epoch,track).mean_replay_rate(index_active_place_fields)];
                data_across_tracks(epoch).track_id=[ data_across_tracks(epoch).track_id track*ones(size(remapping(epoch,track).place_fields_BAYESIAN(index_active_place_fields)))];
            end
            if ~isempty(index_NOTactive_place_fields_that_replay)
                data_across_tracks(epoch).place_field_peak=[data_across_tracks(epoch).place_field_peak remapping(epoch,track).place_fields_BAYESIAN(index_NOTactive_place_fields_that_replay)];
                data_across_tracks(epoch).mean_replay_rate=[data_across_tracks(epoch).mean_replay_rate; remapping(epoch,track).mean_replay_rate(index_NOTactive_place_fields_that_replay)];
                data_across_tracks(epoch).track_id=[ data_across_tracks(epoch).track_id -track*ones(size(remapping(epoch,track).place_fields_BAYESIAN(index_NOTactive_place_fields_that_replay)))];  %id cells without place fields with a negative track value
            end
        end
        
        
    end
    
    
   
    
    figure(current_figure_handle)
    for epoch=1:size(remapping,1)
        subplot(size(remapping,1),1,epoch)
        hold on
        
        all_data.place_field=[];
        all_data.replay=[];
        for track_pair=1:size(remapping,2)
            index=find(data_across_tracks(epoch).track_id==track_pair);
            plot(data_across_tracks(epoch).place_field_peak(index)+1e-1,data_across_tracks(epoch).mean_replay_rate(index)+1e-1,parameters.plot_color_symbol{track_pair});  %add 1e-3 for loglog plot
            
            all_data.place_field=[all_data.place_field data_across_tracks(epoch).place_field_peak(index)];
            all_data.replay=[all_data.replay; data_across_tracks(epoch).mean_replay_rate(index)];
            
            index2=find(data_across_tracks(epoch).track_id==-track_pair);
            plot(data_across_tracks(epoch).place_field_peak(index2)+1,data_across_tracks(epoch).mean_replay_rate(index2)+1e-1,'ks');  %add 1e-3 for loglog plot
            
            
           
        end
        index_non_NaNs=find(~isnan(all_data.place_field));
        if ~isempty(index_non_NaNs)
        lm = fitlm(all_data.place_field(index_non_NaNs)', all_data.replay(index_non_NaNs),'linear')
        [p,F,d] = coefTest(lm);
        linear_fit_P_values_RATE_VS_PEAK_PF=p
        end
        plot([1 10],[.5 5],'k--');
        plot([1 10],[1 10],'k--');
        ax = gca;
        ax.XScale = 'log';
        ax.YScale = 'log';
    end
    
    subplot(2,1,1)
    title('sleep- PRE');
     ylabel('firing rate during replay')
            xlabel('peak response of place field')
    subplot(2,1,2)
    title('sleep- POST');
     ylabel('firing rate during replay')
            xlabel('peak response of place field')
end

end