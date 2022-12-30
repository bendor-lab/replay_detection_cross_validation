function [ROC,log_odd,AUC_combined,AUC_session,Mean_zdata,SE_zdata,pval_average,pval_session] = plot_log_odd_new(ROC,log_odd,timebin,option,figure_option)


% colour_line={'g','k','b','r','m'};
% colour_symbol={'go','ko','bo','ro','mo'};


states=[-1 0 1 2]; % PRE, POST and RUN
% no_of_sessions = 5;
no_of_sessions = max(log_odd.experiment);

%% Replay Index for each replay (sorted by session and state)
session_index = [];
for session = 1:no_of_sessions
    for epoch = 1:length(states) % sort track id and behavioral states
        if isempty(intersect(find(log_odd.experiment == session),find(log_odd.behavioural_state == states(epoch))))
            session_index{session}{epoch} = [];
        else
            session_index{session}{epoch} = intersect(find(log_odd.experiment == session),find(log_odd.behavioural_state == states(epoch)));
        end
    end
end

track_id = log_odd.track;

index = [];% Sort replay id according to session and behavioral states
for session = 1:no_of_sessions
    for k=1:length(states) % sort track id and behavioral states
        state_index = session_index{session}{k};
        
        if isempty(state_index)
            index{session}{1}{k} = [];
            index{session}{2}{k} = [];
        else
            
            if k == 1 % PRE
                index{session}{1}{1} = intersect(state_index,find(track_id==1));
                index{session}{2}{1} = intersect(state_index,find(track_id==2));
            elseif k == 2 % POST
                index{session}{1}{3} = intersect(state_index,find(track_id==1));
                index{session}{2}{3} = intersect(state_index,find(track_id==2));
            elseif k == 3 % RUN Track 1
                index{session}{1}{2} = intersect(state_index,find(track_id==1)); % Track 1 replay on Track 1 (local replay)
            elseif k == 4 % RUN Track 2
                index{session}{2}{2} = intersect(state_index,find(track_id==2)); % Track 2 replay on Track 2 (local replay)
            end
        end
    end
end


% % Resampling so that no of T1 replay matches no of T2 replay 
% for session = 1:no_of_sessions
%     for epoch = 1:3
%         difference = abs(length(index{session}{1}{epoch}) - length(index{session}{2}{epoch})); % find the difference between two tracks
%         
%         if length(index{session}{1}{epoch}) > length(index{session}{2}{epoch})
%             if isempty(index{session}{2}{epoch})
%                 index{session}{1}{epoch} = [];
%             else
%                 index{session}{2}{epoch} = [index{session}{2}{epoch} datasample(index{session}{2}{epoch},difference)]
%             end
%             
%         elseif length(index{session}{1}{epoch}) < length(index{session}{2}{epoch})
%             if isempty(index{session}{1}{epoch}) % If track 1 is empty then track 2 also empty (happens during PRE)
%                 index{session}{2}{epoch} = [];
%             else
%                 index{session}{1}{epoch} = [index{session}{1}{epoch} datasample(index{session}{1}{epoch},difference)]
%             end
%         end
%     end
% end


%% Track ID shuffle
epoch_track_index = [];

for k=1:length(states) % sort track id and behavioral states (pooled from 5 sessions)
    state_index = find(log_odd.behavioural_state==states(k));
    
    if k == 1 % PRE
        epoch_track_index{1}{1} = intersect(state_index,find(track_id==1));
        epoch_track_index{2}{1} = intersect(state_index,find(track_id==2));
    elseif k == 2 % POST
        epoch_track_index{1}{3} = intersect(state_index,find(track_id==1));
        epoch_track_index{2}{3} = intersect(state_index,find(track_id==2));
    elseif k == 3 % RUN Track 1
        epoch_track_index{1}{2} = intersect(state_index,find(track_id==1)); % Track 1 replay on Track 1 considered RUN Replay
%         tempt = intersect(state_index,find(track_id==2)) % Track 2 replay on Track 1 (Preplay)
%         epoch_track_index{2}{1} = [epoch_track_index{2}{1} tempt]; % Put it into PRE              
    elseif k == 4 % RUN Track 2
%         epoch_track_index{1}{3} = intersect(state_index,find(track_id==1));
        epoch_track_index{2}{2} = intersect(state_index,find(track_id==2));        
    end
end

% Resampling so that no of T1 replay matches no of T2 replay 
% for epoch=1:3 % sort track id and behavioral states
%     difference = abs(length(epoch_track_index{1}{epoch}) - length(epoch_track_index{2}{epoch})); % find the difference between two tracks
%     combined = [epoch_track_index{1}{epoch},epoch_track_index{2}{epoch}];
%     
%     if length(epoch_track_index{1}{epoch}) > length(epoch_track_index{2}{epoch})
%         combined = [combined datasample(epoch_track_index{2}{epoch},difference)]; % combine and equate the sample size
%         balanced_index{1}{epoch} = combined(1:length(combined)/2); % divide half into one
%         balanced_index{2}{epoch} = combined(1+length(combined)/2:end); % divide another half into another one
%     elseif length(epoch_track_index{1}{epoch}) < length(epoch_track_index{2}{epoch})
%         if ~isempty(epoch_track_index{1}{epoch}) % If track 1 is empty then track 2 also empty
%             combined = [datasample(epoch_track_index{1}{epoch},difference) combined];
%             balanced_index{1}{epoch} = combined(1:length(combined)/2);
%             balanced_index{2}{epoch} = combined(1+length(combined)/2:end);
%         end
%     end
% end


% 
%% Track ID Shuffled index

shuffled_index = []; % 1000 track Id shuffles across sessions
for shuffle = 1:1000
    for epoch=1:3 % sort track id and behavioral states
        combined = [epoch_track_index{1}{epoch},epoch_track_index{2}{epoch}];
        combined = combined(randperm(length(combined)));
        shuffled_index{shuffle}{1}{epoch} = combined(1:length(epoch_track_index{1}{epoch}));
        shuffled_index{shuffle}{2}{epoch} = combined(1+length(epoch_track_index{1}{epoch}):end);
    end
end

% Resampling so that no of T1 replay matches no of T2 replay 
% shuffled_index = []; % 1000 track Id shuffles across sessions
% for shuffle = 1:1000
%     for epoch=1:3 % sort track id and behavioral states
%         difference = abs(length(epoch_track_index{1}{epoch}) - length(epoch_track_index{2}{epoch})); % find the difference between two tracks
%         combined = [epoch_track_index{1}{epoch},epoch_track_index{2}{epoch}];
%         
%         if length(epoch_track_index{1}{epoch}) > length(epoch_track_index{2}{epoch})
%             combined = [combined datasample(epoch_track_index{2}{epoch},difference)]; % combine and equate the sample size
%             combined = combined(randperm(length(combined))); % track id shuffle
%             shuffled_index{shuffle}{1}{epoch} = combined(1:length(combined)/2); % divide half into one
%             shuffled_index{shuffle}{2}{epoch} = combined(1+length(combined)/2:end); % divide another half into another one
%         elseif length(epoch_track_index{1}{epoch}) < length(epoch_track_index{2}{epoch})
%             if ~isempty(epoch_track_index{1}{epoch}) % If track 1 is empty then track 2 also empty (happens during PRE)
%                 combined = [combined datasample(epoch_track_index{1}{epoch},difference)];
%                 combined = combined(randperm(length(combined)));
%                 shuffled_index{shuffle}{1}{epoch} = combined(1:length(combined)/2);
%                 shuffled_index{shuffle}{2}{epoch} = combined(1+length(combined)/2:end);
%             end
%         end
%     end
% end
% 

% Load data
if strcmp(option,'original')
    if timebin == 1
        
        data = log_odd.one_bin_zscored.original;
%         data(isnan(data)) = 0;
        
    elseif timebin == 0.02
        
        data = log_odd.normal_zscored.original;
    end
    
elseif strcmp(option,'rate fixed')
    if timebin == 1     
        data = log_odd.one_bin_zscored.rate_fixed;
%         data(isnan(data)) = 0;
        
    elseif timebin == 0.02        
        data = log_odd.normal_zscored.rate_fixed;
    end
    
% elseif strcmp(option,'rate fixed median spike count')
%     
%     data = log_odd.one_bin_zscored.rate_fixed_median;
elseif strcmp(option,'normal circular shift')
    data = log_odd.normal_zscored.rate_fixed_circular_shift;
    
elseif strcmp(option,'rate fixed circular shift')
    data = log_odd.normal_zscored.normal_circular_shift;
    
elseif strcmp(option,'rate fixed shuffled spike train')
    data = log_odd.normal_zscored.rate_fixed_shuffled_spike_train;
elseif strcmp(option,'rate fixed global remapping')
    data = log_odd.normal_zscored.rate_fixed_global_remapping;
elseif strcmp(option,'rate remapped one bin') 
    data = log_odd.one_bin_zscored.rate_remapping;      
end


for epoch = 1:3 % PRE, POST, RUN
    
    for nshuffle = 1:1000
        shuffle_average{1}{epoch}(nshuffle) = mean(data(shuffled_index{nshuffle}{1}{epoch})) - mean(data(shuffled_index{nshuffle}{2}{epoch}));
    end
    
    for session= 1:no_of_sessions
        if isempty(data(index{session}{1}{epoch})) % If Preplay is missing T1 replay (ignore it)
            zdata{session}{1}{epoch} = [];
            zdata{session}{2}{epoch} = [];
        else
            
            zdata{session}{1}{epoch} = data(index{session}{1}{epoch});
            zdata{session}{2}{epoch} = data(index{session}{2}{epoch});
            
            % Difference between mean T1 log odds - T2 log odd
            Mean_zdata{1}{session}(1,epoch) = mean(zdata{session}{1}{epoch}) - mean(zdata{session}{2}{epoch});
            % I am not sure how to calculate the SE for the difference between
            % two distribution?
            
            for nshuffle = 1:1000
                tempt = [data(index{session}{1}{epoch}) data(index{session}{2}{epoch})];
                % shuffle replay event track label within the session and
                % within the behavioral epoch
                tempt = tempt(randperm(length(tempt)));
                % Calculate mean T1 (but shuffled) log odd - T2 (but shuffled) log odd
                shuffle_session{session}{epoch}(nshuffle) = mean(tempt(1:length(index{session}{1}{epoch}))) - mean(tempt(1+length(index{session}{1}{epoch}):end));
            end
            
            % Get a zscore and p value based on the distribution of the
            % event track label shuffle for each session.
            tempt = [];
            tempt = (Mean_zdata{1}{session}(1,epoch) - mean(shuffle_session{session}{epoch}))/std(shuffle_session{session}{epoch});
            
            pval_session{1}{session}{epoch} = normcdf(tempt, 0, 1,'upper'); % One tailed test

        end       
    end
    
    tempt1 = [];
    for session = 1:no_of_sessions
        tempt1(session) = Mean_zdata{1}{session}(1,epoch);
        shuffle_average{2}{epoch}(session) = mean(shuffle_session{session}{epoch});
    end

  
%     zdata_difference(:,epoch) = datasample(data(epoch_track_index{1}{epoch}),1000) - datasample(data(epoch_track_index{2}{epoch}),100);

   
    Mean_zdata{2}(epoch) = mean(tempt1);
    SE_zdata(epoch) = std(tempt1)/sqrt(length(tempt1));
% 
%     Mean_zdata{2}(epoch) = mean(data(epoch_track_index{1}{epoch})) - mean(data(epoch_track_index{2}{epoch}));
%     SE_zdata(epoch) = std(tempt1)/sqrt(length(tempt1));

    % Get a zscore and p value based on the comparision between the real mean value and the distribution of the
    % event track label shuffle for each session.
    tempt = [];
    tempt = (Mean_zdata{2}(epoch) - mean(shuffle_average{1}{epoch}))/std(shuffle_average{1}{epoch});
    

    pval_average{1}{epoch} = normcdf(tempt, 0, 1,'upper'); % one tailed
    pval_average{2}{epoch} = ranksum(tempt1,shuffle_average{2}{epoch},'tail','right');

end



if strcmp(figure_option,'Y')
%     state_label = {'PRE T1-T2 log odd difference','RUN T1-T2 log odd difference','POST T1-T2 log odd difference'};
    colour_line= {'k','b','r','g'};
%     colour_line2 = {'k--','b--','r--','g--'};
    colour_symbol={'ko','bo','ro','go'};
%     colour_symbol2= {'kx','bx','rx','gx'};
    colour_fill = {[0 0 0,0.3],[0, 0, 1, 0.3],[1,0,0,0.3]};
    
    nfig = 1;
    fig = figure(nfig)
    fig.Position = [680 300 340 379];
    
    
    for epoch = 1:3
        tempt1 = [];
        for session= 1:no_of_sessions
            if ~isempty(zdata{session}{1}{epoch})
                %             scatter(mean(zdata{session}{1}{epoch}),epoch,colour_symbol{epoch});
                scatter(epoch,Mean_zdata{1}{session}(1,epoch),20,colour_symbol{epoch});
                %             tempt1 = [tempt1 mean(zdata{session}{1}{epoch})];
                hold on
                %             scatter(mean(zdata{session}{2}{epoch}),epoch,colour_symbol2{epoch});
                %             tempt2 = [tempt2 mean(zdata{session}{2}{epoch})];
                %         boxplot(epoch,zdata{session}{2}{epoch})
            end
        end
        tempt1 = Mean_zdata{2}(1,epoch);
        %     tempt2 = mean(data(epoch_track_index{2}{epoch}));
        
        h(epoch) = plot([epoch*ones-0.25 epoch+0.25],[mean(tempt1) mean(tempt1)],colour_line{epoch},'LineWidth',2);
        %     h2(epoch) = plot([mean(tempt2) mean(tempt2)],[epoch*ones-0.25 epoch+0.25],colour_line2{epoch},'DisplayName',state_label2{epoch});
        %     rectangle('Position',[mean(shuffle_average{epoch})-std(shuffle_average{epoch}),epoch-0.1,std(shuffle_average{epoch})*2,0.2],'EdgeColor',colour_line{epoch});
        
%         SEM_shuffle = std(shuffle_average{2}{epoch})/sqrt(length(shuffle_average{2}{epoch}));   
        SEM_shuffle = std(shuffle_average{2}{epoch}); 
        rectangle('Position',[epoch-0.1,mean(shuffle_average{2}{epoch})-SEM_shuffle,0.2,SEM_shuffle*2],...
            'EdgeColor',colour_fill{epoch},'FaceColor', colour_fill{epoch});
        %     rectangle('Position',[mean(shuffle_average{2}{epoch})-std(shuffle_average{2}{epoch}),epoch-0.1,std(shuffle_average{2}{epoch})*2,0.2],'EdgeColor',colour_line{epoch},'LineStyle',':');
        if pval_average{2}{epoch} < 0.05
%             scatter([2.05:1:5.05],AUC(epoch,2:end)+0.03,color_fill{epoch},'*')
            scatter([epoch+0.2],[mean(tempt1)+0.1],colour_line{epoch},'*');
        end
    end
    
%     text1 = sprintf('Mean difference = %.2f (p = %.2e)',Mean_zdata{2}(1,1),pval_average{1}{1})
%     text2 = sprintf('Mean difference = %.2f (p = %.2e)',Mean_zdata{2}(1,2),pval_average{1}{2})
%     text3 = sprintf('Mean difference = %.2f (p = %.2e)',Mean_zdata{2}(1,3),pval_average{1}{3})
%     lgd = legend(h(1:3),{text1,text2,text3},'Location','northoutside')
%     lgd.FontSize = 12;
% %     lgd.Position = [0.1780 0.9288 0.3823 0.0750];
%     legend boxoff
    
    ax = gca;
    set(ax,'LineWidth',1.5)
    ax.YAxis.TickDirection =  'out';       %repeat for XAxis
    ax.YAxis.TickLength =  [.005 1];       %repeat for XAxis
    ax.XAxis.TickDirection =  'out';       %repeat for XAxis
    ax.XAxis.TickLength =  [.005 1];       %repeat for XAxis
    ax.FontSize = 12;

    daspect([1 0.7 1])
%     legend([h1],'location', 'northeast');
    % legend([h1,h2],'location', 'northeast');
    xticks([1,2,3])
    xticklabels({'PRE','RUN','POST'}) 
    ylim( [-1 3.5])
    xlim([0.5 3.5])
    ylabel('z-scored log odds track difference')
    box off
    
    if ~isfolder({'log_odd_figure'})
        mkdir('log_odd_figure')
    end
    
    cd log_odd_figure
%     set(gcf,'PaperSize',[20 10]); %set the paper size to what you want
%     print(gcf,'filename','-dpdf') % then print it
    filename = sprintf('log odds difference %.2i %s.fig',timebin,option)
    saveas(gcf,filename)

    filename = sprintf('log odds difference %.2i %s.pdf',timebin,option)
    saveas(gcf,filename)

    filename = sprintf('log odds difference %.2i %s.svg',timebin,option)
    saveas(gcf,filename)
    cd ..
%     
%     if ~isempty(posbin)
%         if strcmp(option,'original')
%             title({sprintf('Log Odd T1-T2 Difference - %s (%i time bin, 1 positon bins per track)',place_cell_type,timebin)})
%         else
%             title({sprintf('Log Odd T1-T2 Difference (%s) - %s (%i time bin, 1 positon bins per track)',option,place_cell_type,timebin)})
%         end
%     else
%         if strcmp(option,'original')
%             title({sprintf('Log Odd T1-T2 Difference - %s (%i time bin, 20 positon bins per track)',place_cell_type,timebin)})
%         else
%             title({sprintf('Log Odd T1-T2 Difference (%s) - %s (%i time bin, 20 positon bins per track)',option,place_cell_type,timebin)})
%         end
%     end
end


% 
% % Resampling so that no of T1 replay matches no of T2 replay 
% for session = 1:no_of_sessions
%     for epoch = 1:3
%         difference = abs(length(index{session}{1}{epoch}) - length(index{session}{2}{epoch})); % find the difference between two tracks
%         
%         if length(index{session}{1}{epoch}) > length(index{session}{2}{epoch})
%             if isempty(index{session}{2}{epoch})
%                 index{session}{1}{epoch} = [];
%             else
%                 index{session}{2}{epoch} = [index{session}{2}{epoch} datasample(index{session}{2}{epoch},difference)]
%             end
%             
%         elseif length(index{session}{1}{epoch}) < length(index{session}{2}{epoch})
%             if isempty(index{session}{1}{epoch}) % If track 1 is empty then track 2 also empty (happens during PRE)
%                 index{session}{2}{epoch} = [];
%             else
%                 index{session}{1}{epoch} = [index{session}{1}{epoch} datasample(index{session}{1}{epoch},difference)]
%             end
%         end
%     end
% end

% 
% %%
% % ROC analysis Per Session
% %
% %%
% 
% for session = 1:no_of_sessions
%     for epoch = 1:3
%         session_epoch_index{session}{epoch} = [index{session}{1}{epoch} index{session}{2}{epoch}];
%     end
% end
% 
% index = [];
% index = session_epoch_index;
% 
% track_id = zeros(1,length(log_odd.track));
% track_id(find(log_odd.track == 1)) = 1;
% 
% 
% AUC_session = [];
% dprime_session = [];
% FalsePR = [];
% TruePR = [];
% % assume track 1 is YES/1 and Track 2 is NO/0
% 
% ppf = 6;
% 
% for session = 1:no_of_sessions
%     if strcmp(figure_option,'Y')
%         nfigure = nfig + ceil(session/ppf);
%         fig = figure(nfigure);
%         fig.Position = [100 -200 1400 1000];
%         subplot(2,3,session-(nfigure-nfig-1)*ppf);
%     end
% 
%     for epoch = 1:3
%         parfor boots = 1:1000 % Bootstrapping by resampling 1000 times with replacement
%             % Get True and false positive rate and area under the curve
%             resample_index = datasample(index{session}{epoch},length(index{session}{epoch}));
%             [X,Y,T,A,OPTROCPT] = perfcurve(track_id(resample_index),data(resample_index),1,'NBoot',1,'XVals',[0:0.05:1]);
%             
%             AUC(boots) = A(:,1);
%             TPR(:,boots) = Y(:,1);
%         end  
%         
%         AUC_session{session}(epoch,:) = AUC;
%         TruePR{epoch} = TPR;
%         X = [0:0.05:1];
%         FalsePR{epoch}(:) = X;
%         
%         %         dprime_session{session}(epoch,:) = sqrt(2)*norminv(A);
%         
%         if strcmp(figure_option,'Y')
%             if strcmp(option,'rate fixed') && timebin == 1
%                 tempt = [];
%                 tempt = repmat(linspace(TPR(1,1),TPR(end,1),21)',1,1000); % Perfcurve fails to give right True PR when data is all 0
%                 TPR = [];
%                 TPR = tempt;
%             end
%             
%             if timebin == 0.02
%                 if strcmp(option,'original')
%                     ROC(epoch).normal.original(session).session = session;
%                     ROC(epoch).normal.original(session).FalsePR = X;
%                     ROC(epoch).normal.original(session).TruePR = TPR;
%                     ROC(epoch).normal.original(session).AUC = AUC;
%                 elseif strcmp(option,'rate fixed')
%                     ROC(epoch).normal.rate_fixed(session).session = session;
%                     ROC(epoch).normal.rate_fixed(session).FalsePR = X;
%                     ROC(epoch).normal.rate_fixed(session).TruePR = TPR;
%                     ROC(epoch).normal.rate_fixed(session).AUC = AUC;
%                 elseif strcmp(option,'rate fixed global remapping')
%                     ROC(epoch).normal.rate_fixed_global_remapped(session).session = session;
%                     ROC(epoch).normal.rate_fixed_global_remapped(session).FalsePR = X;
%                     ROC(epoch).normal.rate_fixed_global_remapped(session).TruePR = TPR;
%                     ROC(epoch).normal.rate_fixed_global_remapped(session).AUC = AUC;
%                 end
%                 
%             elseif timebin == 1
%                 if strcmp(option,'original')
%                     ROC(epoch).one_bin.original(session).session = session;
%                     ROC(epoch).one_bin.original(session).FalsePR = X;
%                     ROC(epoch).one_bin.original(session).TruePR = TPR;
%                     ROC(epoch).one_bin.original(session).AUC = AUC;
%                 elseif strcmp(option,'rate fixed')
%                     ROC(epoch).one_bin.rate_fixed(session).session = session;
%                     ROC(epoch).one_bin.rate_fixed(session).FalsePR = X;
%                     % Because Perfcurve function fails to give right y when all log odds is zero.
%                     ROC(epoch).one_bin.rate_fixed(session).TruePR = repmat(linspace(TPR(1,1),TPR(end,1),21)',1,1000);
%                     ROC(epoch).one_bin.rate_fixed(session).AUC = AUC;
%                 elseif strcmp(option,'rate remapped one bin')
%                     ROC(epoch).one_bin.rate_remapping(session).session = session;
%                     ROC(epoch).one_bin.rate_remapping(session).FalsePR = X;
%                     ROC(epoch).one_bin.rate_remapping(session).TruePR = TPR;
%                     ROC(epoch).one_bin.rate_remapping(session).AUC = AUC;
%                 end
%             end
%             
%             % Calculate SEM from 95% confidence interval
%             p(epoch)= plot(X,mean(TPR,2),colour_line{epoch});
%             
%             for i = 1:size(TPR,1)
%                 Y_SE(i) = std(TPR(i,:)); % SD of bootstrapped distribution is Bootstrap SE
%             end
%             
%             hold on
%             % If confidence interval
% %             x2 = [X', fliplr(X')];
% %             inBetween = [Y(:,2)', fliplr(Y(:,3)')];
% %             fill(x2, inBetween, colour_line{epoch},'FaceAlpha','0.25','LineStyle','none');
% %             
%             % If SEM
%             x2 = [X, fliplr(X)];
%             inBetween = [mean(TPR,2)'+Y_SE, fliplr(mean(TPR,2)'-Y_SE)];
%             fill(x2, inBetween, colour_line{epoch},'FaceAlpha','0.25','LineStyle','none');            
%             
%             xlabel('False Positive Rate')
%             ylabel('True positive rate')
%             
%             xlim([0 1])
%             ylim([0 1])
%             plot([0 1],[0 1],'k--');
%             box off
%             
%             AUC = [];
%             TRP = [];
%             FalsePR = [];
%             TruePR= [];
%             TPR = [];
%         end
%     end
%     
%     text1 = sprintf('Mean AUC %.2f',mean(AUC_session{session}(1,:)));
%     text2 = sprintf('Mean AUC %.2f',mean(AUC_session{session}(2,:)));
%     text3 = sprintf('Mean AUC %.2f',mean(AUC_session{session}(3,:)));    
%     lgd = legend(p(1:3),{text1,text2,text3},'Location','southeast');
%     lgd.FontSize = 12;
%     legend boxoff
%     
%     ax = gca;
%     set(ax,'LineWidth',1.5)
%     ax.YAxis.TickDirection =  'out';       %repeat for XAxis
%     ax.YAxis.TickLength =  [.005 1];       %repeat for XAxis
%     ax.XAxis.TickDirection =  'out';       %repeat for XAxis
%     ax.XAxis.TickLength =  [.005 1];       %repeat for XAxis
%     ax.FontSize = 12;
%     axis square
%     
% 
% end
% 
%     cd log_odd_figure
%     figure(2);
%     filename = sprintf('ROC per session %.2i %s 1-6.pdf',timebin,option)
%     set(gcf,'PaperType','A2'); %set the paper size to what you want
%     print(gcf,filename,'-dpdf') % then print it
%     saveas(gcf,filename)
%     
%     figure(3);
%     filename = sprintf('ROC per session %.2i %s 7-10.pdf',timebin,option)
%     set(gcf,'PaperType','A2'); %set the paper size to what you want
%     print(gcf,filename,'-dpdf') % then print it
%     saveas(gcf,filename)    
%     cd ..
% 
% %     str{epoch} = sprintf(['%s mean AUC - %.2f & mean dprime - %.2f'],Legend{epoch},mean(AUC(epoch,:)),mean(dprime(epoch,:)));
% %     text(FalsePR{epoch}{c}(ceil(length(TruePR{epoch}{c})/2)),TruePR{epoch}{c}(ceil(length(TruePR{epoch}{c})/2)),str{epoch});
% % if strcmp(figure_option,'Y')
% %     if ~isempty(posbin)
% %         if strcmp(option,'original')
% %             sgtitle({sprintf('ROC for Classification by Logistic Regression - %s (%i time bin, 1 positon bins per track)',place_cell_type,timebin)})
% %         else
% %             sgtitle({sprintf('ROC for Classification by Logistic Regression(%s) - %s (%i time bin, 1 positon bins per track)',option,place_cell_type,timebin)})
% %         end
% %     else
% %         if strcmp(option,'original')
% %             sgtitle({sprintf('ROC for Classification by Logistic Regression - %s (%i time bin, 20 positon bins per track)',place_cell_type,timebin)})
% %         else
% %             sgtitle({sprintf('ROC for Classification by Logistic Regression(%s) - %s (%i time bin, 20 positon bins per track)',option,place_cell_type,timebin)})
% %         end
% %     end
% % end
% 

%%
% ROC analysis using Session combined
%
%%
% 
% % Resampling so that no of T1 replay matches no of T2 replay 
% for epoch=1:3 % sort track id and behavioral states
%     difference = abs(length(epoch_track_index{1}{epoch}) - length(epoch_track_index{2}{epoch})); % find the difference between two tracks
%     combined = [epoch_track_index{1}{epoch},epoch_track_index{2}{epoch}];
%     
%     if length(epoch_track_index{1}{epoch}) > length(epoch_track_index{2}{epoch})
%         combined = [combined datasample(epoch_track_index{2}{epoch},difference)]; % combine and equate the sample size
%         balanced_index{1}{epoch} = combined(1:length(combined)/2); % divide half into one
%         balanced_index{2}{epoch} = combined(1+length(combined)/2:end); % divide another half into another one
%     elseif length(epoch_track_index{1}{epoch}) < length(epoch_track_index{2}{epoch})
%         if ~isempty(epoch_track_index{1}{epoch}) % If track 1 is empty then track 2 also empty
%             combined = [datasample(epoch_track_index{1}{epoch},difference) combined];
%             balanced_index{1}{epoch} = combined(1:length(combined)/2);
%             balanced_index{2}{epoch} = combined(1+length(combined)/2:end);
%         end
%     end
% end

% 
% epoch_index = [];
% for epoch = 1:3
%     epoch_index{epoch} = [balanced_index{1}{epoch} balanced_index{2}{epoch}];
% end

epoch_index = [];
for k=1:length(states)
    state_index = find(log_odd.behavioural_state==states(k));
    
    if states(k) == 1 % On track 1
        epoch_index{2} = [intersect(state_index,find(log_odd.track==1))]; % Track 1 replay on Track 1 considered RUN Replay
    elseif states(k) == 2 % On track 2
        tempt = intersect(state_index,find(log_odd.track==2));
         epoch_index{2} =  [epoch_index{2} tempt];
    elseif k == 2
        epoch_index{3} = [intersect(state_index,find(log_odd.track==1)) intersect(state_index,find(log_odd.track==2))];
    elseif k == 1
        epoch_index{1} = [intersect(state_index,find(log_odd.track==1)) intersect(state_index,find(log_odd.track==2))];
    end
end

index = epoch_index;
track_id = zeros(1,length(log_odd.track));
track_id(find(log_odd.track == 1)) = 1;


if strcmp(figure_option,'Y')
    nfig = 4;
    fig = figure(nfig)
    fig.Position = [460 290 314 257]
    Legend = {'PRE','RUN','POST'};
end

AUC = [];
dprime = [];
FalsePR = [];
TruePR = [];
AUC_combined = zeros(3,1000);
% track 1 is YES/1 and Track 2 is NO/0


for epoch = 1:3
    tic
    AUC = [];
    TPR = [];
    parfor boots = 1:1000 % Bootstrapping by resampling 1000 times with replacement
        % Get True and false positive rate and area under the curve
        sprintf('Bootstrpping %i',boots)
        resample_index = datasample(index{epoch},length(index{epoch}));
        [X,Y,T,A,OPTROCPT] = perfcurve(track_id(resample_index),data(resample_index),1,'NBoot',1,'XVals',[0:0.05:1]);

        AUC(boots) = A(:,1);
        TPR(:,boots) = Y(:,1);
    end
    toc
    AUC_combined(epoch,:) = AUC;
% end    
%     TruePR{epoch} = TPR;
    X = [0:0.05:1];
%     FalsePR{epoch}(:) = X;

%     dprime(epoch,:) = sqrt(2)*norminv(A);
      
    
    if strcmp(figure_option,'Y')
        if strcmp(option,'rate fixed') && timebin == 1
            tempt = [];
            tempt = repmat(linspace(TPR(1,1),TPR(end,1),21)',1,1000); % Perfcurve fails to give right True PR when data is all 0
            TPR = [];
            TPR = tempt;
        end
        
        % Save ROC curve
        if timebin == 0.02
            if strcmp(option,'original')
                ROC(epoch).normal.original(6).session = 'Combined';
                ROC(epoch).normal.original(6).FalsePR = X;
                ROC(epoch).normal.original(6).TruePR = TPR;
                ROC(epoch).normal.original(6).AUC = AUC;
            elseif strcmp(option,'rate fixed')
                ROC(epoch).normal.rate_fixed(6).session = 'Combined';
                ROC(epoch).normal.rate_fixed(6).FalsePR = X;
                ROC(epoch).normal.rate_fixed(6).TruePR = TPR;
                ROC(epoch).normal.rate_fixed(6).AUC = AUC;
            elseif strcmp(option,'rate fixed global remapping')
                ROC(epoch).normal.rate_fixed_global_remapped(6).session = 'Combined';
                ROC(epoch).normal.rate_fixed_global_remapped(6).FalsePR = X;
                ROC(epoch).normal.rate_fixed_global_remapped(6).TruePR = TPR;
                ROC(epoch).normal.rate_fixed_global_remapped(6).AUC = AUC;
            end
            
        elseif timebin == 1
            if strcmp(option,'original')
                ROC(epoch).one_bin.original(6).session = 'Combined';
                ROC(epoch).one_bin.original(6).FalsePR = X;
                ROC(epoch).one_bin.original(6).TruePR = TPR;
                ROC(epoch).one_bin.original(6).AUC = AUC;
            elseif strcmp(option,'rate fixed')
                ROC(epoch).one_bin.rate_fixed(6).session = 'Combined';
                ROC(epoch).one_bin.rate_fixed(6).FalsePR = X;
                % Because Perfcurve function fails to give right y when all log odds is zero.
                ROC(epoch).one_bin.rate_fixed(6).TruePR = repmat(linspace(TPR(1,1),TPR(end,1),21)',1,1000);
                ROC(epoch).one_bin.rate_fixed(6).AUC = AUC;
            elseif strcmp(option,'rate remapped one bin')
                ROC(epoch).one_bin.rate_remapping(6).session = 'Combined';
                ROC(epoch).one_bin.rate_remapping(6).FalsePR = X;
                ROC(epoch).one_bin.rate_remapping(6).TruePR = TPR;
                ROC(epoch).one_bin.rate_remapping(6).AUC = AUC;
            end
        end
        
        
        p(epoch)= plot(X,mean(TPR,2),colour_line{epoch});
        
        for i = 1:size(TPR,1)
            Y_SE(i) = std(TPR(i,:)); % SD of bootstrapped distribution is Bootstrap SE
        end
        
        hold on
        
        % If confidence interval
        %             x2 = [X', fliplr(X')];
        %             inBetween = [Y(:,2)', fliplr(Y(:,3)')];
        %             fill(x2, inBetween, colour_line{epoch},'FaceAlpha','0.25','LineStyle','none');
        %
        
        % If SEM
        x2 = [X, fliplr(X)];
        inBetween = [mean(TPR,2)'+Y_SE, fliplr(mean(TPR,2)'-Y_SE)];
        fill(x2, inBetween, colour_line{epoch},'FaceAlpha','0.25','LineStyle','none');
        xlabel('False Positive Rate')
        ylabel('True positive rate')
        
        %         str{epoch} = sprintf(['mean AUC - %.3f\nmean dprime - %.3f'],mean(AUC(epoch,:)),mean(dprime(epoch,:)));
        % %         str{epoch} = sprintf(['%s\nmean AUC - %.3f\nmean dprime - %.3f'],Legend{epoch},mean(AUC(epoch,:)),mean(dprime(epoch,:)));
        %         if epoch == 1
        %             text(FalsePR{epoch}{c}(ceil(length(TruePR{epoch}{c})/2)),TruePR{epoch}{c}(ceil(length(TruePR{epoch}{c})/2)),str{epoch},'Color',colour_line{epoch});
        %         elseif epoch == 2
        %             text(FalsePR{epoch}{c}(ceil(length(TruePR{epoch}{c})/2)),[TruePR{epoch}{c}(ceil(length(TruePR{epoch}{c})/2)) - 0.05],str{epoch},'Color',colour_line{epoch});
        %         elseif epoch == 3
        %             text(FalsePR{epoch}{c}(ceil(length(TruePR{epoch}{c})/2)),[TruePR{epoch}{c}(ceil(length(TruePR{epoch}{c})/2))- 0.05],str{epoch},'Color',colour_line{epoch});
        %         end
        %
        
        xlim([0 1])
        ylim([0 1])
        plot([0 1],[0 1],'k--');
        box off
    end


end

text1 = sprintf('AUC %.2f',mean(AUC_combined(1,:)));
text2 = sprintf('AUC %.2f',mean(AUC_combined(2,:)));
text3 = sprintf('AUC %.2f',mean(AUC_combined(3,:)));
lgd = legend(p(1:3),{text1,text2,text3},'Location','southeast');
lgd.FontSize = 12;
legend boxoff

ax = gca;
set(ax,'LineWidth',1.5)
ax.YAxis.TickDirection =  'out';       %repeat for XAxis
ax.YAxis.TickLength =  [.005 1];       %repeat for XAxis
ax.XAxis.TickDirection =  'out';       
ax.XAxis.TickLength =  [.005 1];       
ax.FontSize = 12;

cd log_odd_figure
filename = sprintf('ROC session combined %.2i %s.fig',timebin,option)
saveas(gcf,filename)

filename = sprintf('ROC session combined %.2i %s.pdf',timebin,option)
saveas(gcf,filename)

filename = sprintf('ROC session combined %.2i %s.svg',timebin,option)
saveas(gcf,filename)
cd ..



end
