function plot_PSDs
load PSD_data

if ~exist('./figures/PSD','dir')
    mkdir('figures/PSD');
end

for i=1: length(all_PSD)

    f=figure(i);
    if isfield(all_PSD,'PSD_F_awake')
        ax= subplot(1,2,1);
        title(strcat(['channel: ' num2str(all_PSD(i).CSCchannel) ' awake'])); hold on
        u= plot(all_PSD(i).PSD_F_awake,mean(all_PSD(i).PSD_awake),'k'); hold on
        ax.XScale = 'log';
        xlim([0 max(all_PSD(i).PSD_F_awake)]);
        ylim([0 max(max(mean(all_PSD(i).PSD_awake)),max(mean(all_PSD(i).PSD_sleep)))]);   yl=ylim;    
        area([5 12],[1 1]*yl(2),'FaceColor',[.9 .9 .8],'FaceAlpha',0.8); 
        area([9 17],[1 1]*yl(2),'FaceColor',[.8 .8 .8],'FaceAlpha',0.8);  
        area([40 100],[1 1]*yl(2),'FaceColor',[.9 .8 .9],'FaceAlpha',0.8); 
        area([125 200],[1 1]*yl(2),'FaceColor',[.8 .9 .9],'FaceAlpha',0.8);
        uistack(u,'top');
        legend({'theta','spindle','high gamma','ripple','PSD'},'Location','Southwest');
    
        ax2= subplot(1,2,2);
        title(strcat(['channel: ' num2str(all_PSD(i).CSCchannel) ' sleep'])); hold on
        u2= plot(all_PSD(i).PSD_F_sleep,mean(all_PSD(i).PSD_sleep),'k'); hold on
        ax2.XScale = 'log';   
        xlim([0 max(all_PSD(i).PSD_F_sleep)]);
        ylim([0 max(max(mean(all_PSD(i).PSD_awake)),max(mean(all_PSD(i).PSD_sleep)))]);  yl= ylim;  
        area([5 12],[1 1]*yl(2),'FaceColor',[.9 .9 .8],'FaceAlpha',0.8); 
        area([9 17],[1 1]*yl(2),'FaceColor',[.8 .8 .8],'FaceAlpha',0.8); 
        area([40 100],[1 1]*yl(2),'FaceColor',[.9 .8 .9],'FaceAlpha',0.8);  
        area([125 200],[1 1]*yl(2),'FaceColor',[.8 .9 .9],'FaceAlpha',0.8);
        uistack(u2,'top');
    else
        title(strcat(['channel: ' num2str(all_PSD(i).CSCchannel)]));
        u= plot(all_PSD(i).PSD_F,mean(all_PSD(i).PSD),'k'); hold on
        ax.XScale = 'log';
        xlim([0 300]);
        ylim([0 max(max(mean(all_PSD(i).PSD)),max(mean(all_PSD(i).PSD)))]);   yl=ylim;
        area([5 12],[1 1]*yl(2),'FaceColor',[.9 .9 .8],'FaceAlpha',0.8); 
        area([9 17],[1 1]*yl(2),'FaceColor',[.8 .8 .8],'FaceAlpha',0.8); 
        area([40 100],[1 1]*yl(2),'FaceColor',[.9 .8 .9],'FaceAlpha',0.8);  
        area([125 200],[1 1]*yl(2),'FaceColor',[.8 .9 .9],'FaceAlpha',0.8);
        uistack(u,'top');
        legend({'theta','spindle','high gamma','ripple','PSD'},'Location','Southwest');
    end
    
    saveas(f,['figures/PSD/PSD_CSC' num2str(all_PSD(i).CSCchannel) '.jpg']);
    close(f);
end

end