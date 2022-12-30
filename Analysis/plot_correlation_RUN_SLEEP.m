% Plots correlation between RUN and SLEEP replay dynamics
% INPUTS:
%               'scoring': 'wcorr','spearman' - in the future will add 'linear' (default: 'wcorr')                                    
%               'control': 'fixed_rate','fixed_spike_count','
%                   global_remapped','rate_remapped','rate_detection_control' (default: none)
%               'epochs': list of epochs to include in plot as a cell array of epoch names. need to match availability in
%                   remapping_[...].mat file.  (default: 'PRE' + 'POST')
%               'target_epoch': epoch to compare against
%               'vars': name of fieldname for data plotted on x & y axis, need to match with a field in remapping or remappin_raw              
%               'cdata': cell array of RGB/else vectors for each included cell

function plot_correlation_RUN_SLEEP(varargin)


available_controls= {'fixed_rate','fixed_spike_count','global_remapped','rate_remapped',...
                                'rate_detection_control','rate_intrinsic_bias_control',...
                                'replay_rate_shuffle_control','replay_rate_shuffle_detection'};
available_scoring= {'wcorr','spearman'};

p= inputParser;
addParameter(p,'scoring','wcorr',@(x) ismember(x,available_scoring));
addParameter(p,'epochs',{'PRE','POST'},@iscell);
addParameter(p,'target_epoch',{'RUN'},@iscell);
addParameter(p,'vars',{'track1_median_replay_rate','track2_median_replay_rate'},@(x) ischar(x) || iscell(x));
addParameter(p,'control',[],@(x) ismember(x,available_controls));
addParameter(p,'cdata',[0.5 0.5 0.5],@ismatrix);
addParameter(p,'axis_label',[],@ischar);
parse(p,varargin{:});

filename= ['rate_remapping_analysis_TRACK_PAIRS_' p.Results.scoring];
if isempty(p.Results.control)
    load([filename '.mat']);
else
    switch p.Results.control
        case 'fixed_rate'
            load(['.\CONTROLS\' p.Results.control '\' filename '_FIXED.mat']);
        case 'rate_detection_control'
            load(['.\CONTROLS\rate_remapped\' filename '_RATE.mat']);
        case 'rate_intrinsic_bias_control'
            load(['.\CONTROLS\rate_remapped\' filename '_INTRINSIC_RATE.mat']);
        case 'rate_remapped'
            load(['.\CONTROLS\rate_remapped\' filename '.mat']);
        case 'replay_rate_shuffle_control'
            load(['.\CONTROLS\replay_rate_shuffle\' filename '_REPLAY_RATE.mat']);
        case 'replay_rate_shuffle_detection'
            load(['.\CONTROLS\replay_rate_shuffle\' filename '_REPLAY_RATE_DETECTION.mat']);
        otherwise
            load(['.\CONTROLS\' p.Results.control '\' filename '.mat']);
    end
end

xy_var = p.Results.vars;
if ischar(xy_var)
    xy_var= {xy_var};
end
cdata= p.Results.cdata;
epochs= p.Results.epochs;
[~,epoch_idx]= ismember(epochs,vertcat(remapping.epoch));
target_epoch_idx= find(strcmp([remapping.epoch],p.Results.target_epoch));

f1= figure('units','normalized','Color','w');
f1.Name = [remapping(target_epoch_idx).epoch{:} '-SLEEP track correlation'];
f2= figure('units','normalized','Color','w');
f2.Name = [remapping(target_epoch_idx).epoch{:}  '-SLEEP diff correlation'];
subplot_counter = 1;
subplot_counter2 = 1;

step_size = 2;
edge_val = 15;

for ep =  1 : length(epochs)
    
    if size(remapping_raw(epoch_idx(ep)).new_ID,2) > size(remapping_raw(target_epoch_idx).new_ID,2) %if there are more cells in RUN than in the other field (i.e. PRE/POST)
        [~,idx1,idx2] = intersect(remapping_raw(epoch_idx(ep)).new_ID,remapping_raw(target_epoch_idx).new_ID);
        % Prepare colormap
        min_val = -edge_val; % min(remapping(epoch_idx(ep)).place_field_diff(idx1));
        max_val = edge_val;  % max(remapping(epoch_idx(ep)).place_field_diff(idx1));
        steps = min_val:step_size:max_val;
        steps(end+1) = steps(end) + step_size;
        temp=remapping(epoch_idx(ep)).place_field_diff(idx1);
        temp(temp<min_val) = min_val;
        temp(temp>max_val) = max_val;
        [~,~,bins] = histcounts(temp,steps);
        %[~,~,bins] = histcounts(remapping(epoch_idx(ep)).place_field_diff(idx1),steps);
    else
        [~,idx2,idx1] = intersect(remapping_raw(target_epoch_idx).new_ID,remapping_raw(epoch_idx(ep)).new_ID);
        % Prepare colormap
        min_val = -edge_val; % min(remapping(target_epoch_idx).place_field_diff(idx1));
        max_val = edge_val;  % max(remapping(target_epoch_idx).place_field_diff(idx1));
        steps = min_val:step_size:max_val;
        steps(end+1) = steps(end) + step_size;
        temp=remapping(target_epoch_idx).place_field_diff(idx2);
        temp(temp<min_val) = min_val;
        temp(temp>max_val) = max_val;
        [~,~,bins] = histcounts(temp,steps);
    end
    
    for this_var_pair = 1 : length(xy_var)
        
        if any(~contains(xy_var{this_var_pair},'diff'))
            
            % For each track, compare max instantaneous FR between RUN & POST
            [rho,pval] = corr(remapping_raw(epoch_idx(ep)).(xy_var{this_var_pair})(idx1),remapping_raw(target_epoch_idx).(xy_var{this_var_pair})(idx2));

            
            figure(f1)
            ax(subplot_counter)= subplot(length(epochs),length(find(~contains(xy_var,'diff')==1)),subplot_counter);
            hold(ax(subplot_counter),'on');
            
            scatter(ax(subplot_counter),remapping_raw(epoch_idx(ep)).(xy_var{this_var_pair})(idx1),remapping_raw(target_epoch_idx).(xy_var{this_var_pair})(idx2),...
                20,cdata,'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
            
            if isempty(p.Results.axis_label)
                axlabel_stripped= xy_var{this_var_pair};
                axlabel_stripped(axlabel_stripped == '_')= ' ';
                xlabel([remapping_raw(epoch_idx(ep)).epoch{1,1} ' - ' axlabel_stripped]);
                ylabel([remapping(target_epoch_idx).epoch{:} ' - ' axlabel_stripped]);
            else
                xlabel([remapping_raw(epoch_idx(ep)).epoch{1,1} ' - ' p.Results.axis_label]);
                ylabel([remapping(target_epoch_idx).epoch{:}  ' - ' p.Results.axis_label]);
            end
            
            title(['R = ' num2str(rho) '  P = ' num2str(pval)])
            subplot_counter = subplot_counter+1;
            
        elseif any(contains(xy_var{this_var_pair},'diff'))
            
            % For each track, compare max instantaneous FR between RUN & POST
              %[rho,pval] = corr(remapping(epoch_idx(ep)).(xy_var{this_var_pair})(idx1),remapping(target_epoch_idx).(xy_var{this_var_pair})(idx2));
            lm = fitlm(remapping(epoch_idx(ep)).(xy_var{this_var_pair})(idx1),remapping(target_epoch_idx).(xy_var{this_var_pair})(idx2),'linear');
            [pval_lm,F_lm,~] = coefTest(lm);
            x=[min(remapping(epoch_idx(ep)).(xy_var{this_var_pair})(idx1)) max(remapping(epoch_idx(ep)).(xy_var{this_var_pair})(idx1))];
            b=lm.Coefficients.Estimate';
            slope_intcpt_lm = fliplr(b); 
            disp(['Pval=' num2str(pval_lm) ' ;F_lm=' num2str(F_lm) ' ;slope_intcpt_lm=' num2str(slope_intcpt_lm(1)) ' ;B=' num2str(b(1))])
            % To get coefficients
              %coefs = polyfit(remapping(target_epoch_idx).(xy_var{this_var_pair})(idx2),remapping(epoch_idx(ep)).(xy_var{this_var_pair})(idx1), 1);
              %coefs(1)
            
            % Colormap
            %colmap = jet(length(steps));
            
            % Custom made colormap
            colmap = zeros(length(steps),3);
            colmap(:,1) = linspace(0,1,length(steps)); %0,1
            colmap(:,2) = linspace(0,1,length(steps));
            colmap(:,3) = linspace(1,0,length(steps)); %1,0
            cdata_diff = zeros(length(bins),3);
            for k = 1 : length(steps)
                cdata_diff(bins == k,:) = repmat(colmap(k,:),length(find(bins == k)==1),1);
            end
            
            
            figure(f2)
            ax2(subplot_counter2)= subplot(length(epochs),length(find(contains(xy_var,'diff')==1)),subplot_counter2);
            hold(ax2(subplot_counter2),'on');
            
            scatter(ax2(subplot_counter2),remapping(epoch_idx(ep)).(xy_var{this_var_pair})(idx1),remapping(target_epoch_idx).(xy_var{this_var_pair})(idx2),...
                20,cdata_diff,'filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1);
            caxis([min_val max_val])
            h(subplot_counter2) = colorbar;  
            colormap(colmap)
                        
            plot(x,polyval(fliplr(b),x),'k--');
            
            if isempty(p.Results.axis_label)
                axlabel_stripped= xy_var{this_var_pair};
                axlabel_stripped(axlabel_stripped == '_')= ' ';
                %xlabel([remapping_raw(epoch_idx(ep)).epoch{1,1} ' - ' axlabel_stripped]);
               % ylabel([remapping(target_epoch_idx).epoch{:} ' - ' axlabel_stripped]);
                xlabel(axlabel_stripped);
                ylabel(axlabel_stripped);
            else
                %xlabel([remapping_raw(epoch_idx(ep)).epoch{1,1} ' - ' p.Results.axis_label]);
                %xlabel(p.Results.axis_label);
                %ylabel([remapping(target_epoch_idx).epoch{:} ' - ' p.Results.axis_label]);
                xlabel(p.Results.axis_label);
                ylabel(p.Results.axis_label);
            end
            
            
            %title(['R = ' num2str(rho) '  P = ' num2str(pval)])
            %title(['F = ' num2str(F_lm) '  P = ' num2str(pval_lm)])
            text(ax2(subplot_counter2),.7,0.1,['p = ' num2str(pval_lm)],'Units','Normalized','FontSize',12,'FontName','Arial');
            
            subplot_counter2 = subplot_counter2+1;
            
        end
    end
end

run_format_settings(f2,'match_ax') %for this plot, change formatting_settings NumYTicks to 5

% if any(ismember(xy_var,'mean_max_FR_replay_diff'))
%     allax = flipud(findall(f2,'type','axes'));
%     ib = find(ismember(xy_var,'mean_max_FR_replay_diff'))*[1:length(epochs)];
%     lnk_ay = allax(ib);
%     for j = 1 : length(lnk_ay)
%         yt = cellfun(@str2num,yticklabels).*sum(gausswin(0.25/(1/1000)));
%         lnk_ay(j).YTickLabels = round(yt,2);
%         xt = cellfun(@str2num,xticklabels).*sum(gausswin(0.25/(1/1000)));
%         lnk_ay(j).XTickLabels = round(xt,2);
%     end
% end


if ~isempty(f2)
    for j = 1 : length(ax2)
        ax2(j).XLabel.String = [remapping_raw(epoch_idx(j)).epoch{1,1} ' - ' ax2(j).XLabel.String];
        ax2(j).YLabel.String = [remapping(target_epoch_idx).epoch{:} ' - ' ax2(j).YLabel.String];
    end
    if ~isempty(h)
        for j = 1 : length(h)
            h(j).TickLabels{1} = '< 15';
            h(j).TickLabels{end} = '> 15';
        end
    end
end
end