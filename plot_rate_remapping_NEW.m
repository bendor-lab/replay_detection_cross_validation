function [pval,F,slope_intcpt_lm]= plot_rate_remapping_NEW(varargin)
%%% now as a name-var pair input method
% INPUTS:
%               'scoring': 'wcorr','spearman' - in the future will add 'linear' (default: 'wcorr')                                    
%               'method': 'TRACK_PAIRS' or 'ONE_TRACK' (default: 'TRACK_PAIRS')        
%               'control': 'fixed_rate','fixed_spike_count','
%                   global_remapped','rate_remapped','rate_detection_control' (default: none)
%               'epochs': list of epochs to include in plot as a cell array of epoch names. need to match availability in
%                   remapping_[...].mat file.  (default: 'PRE' + 'POST')
%               'common_cells_toggle': use only cells common to epochs (default: 0)
%               'x_var': name of fieldname for data plotted on x axis
%               'y_var': name of fieldname for data plotted on y axis
%               'subset': fieldname from subsets_of_cells.mat
%               'cdata': cell array of RGB/else vectors for each included cell
%               'x_label': string
%               'y_label': string
%               'colormap': nx3 colormap
%               'colorbar': 0 or 1

load('X:\BendorLab\Drobo\Neural and Behavioural Data\Rate remapping\Data\Tables\subsets_of_cells.mat');
available_controls= {'fixed_rate','fixed_spike_count','global_remapped','rate_remapped',...
                                'rate_detection_control','rate_intrinsic_bias_control',...
                                'replay_rate_shuffle_control','replay_rate_shuffle_detection'};
available_scoring= {'wcorr','spearman'};
available_methods= {'TRACK_PAIRS','ONE_TRACK'};
available_subsets= unique(subset_of_cells.subset)';

p= inputParser;
addParameter(p,'use_mat',[],@isstruct);
addParameter(p,'scoring','wcorr',@(x) ismember(x,available_scoring));
addParameter(p,'control',[],@(x) ismember(x,available_controls));
addParameter(p,'method','TRACK_PAIRS',@(x) ismember(x,available_methods));
addParameter(p,'epochs',{'PRE','POST'},@iscell);
addParameter(p,'x_var',{'place_field_diff'},@(x) ischar(x) || iscell(x));
addParameter(p,'y_var',{'mean_max_FR_replay_diff'},@(x) ischar(x) || iscell(x));
addParameter(p,'common_cells_toggle',0,@ismatrix);
addParameter(p,'subset',[],@(x) any(contains(x,available_subsets)));
addParameter(p,'cdata',[0.5 0.5 0.5],@ismatrix); %[0.2 0.4 0.8]
addParameter(p,'x_label','',@(x) ischar(x) || iscell(x)); % had to change to not overlap with function name
addParameter(p,'y_label','',@(x) ischar(x) || iscell(x));
addParameter(p,'colormap',[],@ismatrix);
addParameter(p,'colorbar',0,@ismatrix);
addParameter(p,'session',0,@ismatrix); % default is 0 for all sessions, otherwise number from 1 to 5
parse(p,varargin{:});

if isempty(p.Results.use_mat)
    filename= ['rate_remapping_analysis_' p.Results.method '_' p.Results.scoring];
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
else
    remapping= p.Results.use_mat;
end

if ~isempty(p.Results.subset)
    if iscell(p.Results.subset) & length(p.Results.subset) > 1
        subset = intersect(subset_of_cells.cell_IDs{ismember(subset_of_cells.subset,p.Results.subset)}); %multintersect?
    else
        
        subset= subset_of_cells.cell_IDs{strcmp(subset_of_cells.subset,p.Results.subset)};
    end
else
    subset= [];
end


% cdata= p.Results.cdata;
experiment_idx = p.Results.session;
epochs= p.Results.epochs;
[~,epoch_idx]= ismember(epochs,vertcat(remapping.epoch));
x_var= p.Results.x_var;
y_var= p.Results.y_var;
if ischar(x_var)
    x_var= {x_var};
end
if ischar(y_var)
    x_var= {y_var};
end


%% create figure
% one row per epoch, one column per [x_var, y_var] pair

for track_pair = 1:size(remapping,2)
    f1= figure('units','normalized','Color','w');
    % f1= figure('Position',[662.3333333333333,591,830,646.6666666666665],'Color',[1 1 1]);
%     current_figure_handle=gcf;
    parameters = list_of_parameters;

    foldername = strsplit(pwd,'\');
    f1.Name= strcat(p.Results.scoring,'_',foldername{end});

    subplot_counter=1;
    % make sure repetitions are not discarded
    if length(epoch_idx)>1
        common_epoch_cells= multintersect(remapping(epoch_idx,track_pair).new_ID);
    else
        common_epoch_cells= remapping(epoch_idx,track_pair).new_ID;
    end

    for this_epoch=1:length(epochs)
        for this_var_pair=1:length(x_var)
        
            % check the var names are fields from loaded filename
            if ~isfield(remapping,x_var{this_var_pair}) ||  ~isfield(remapping,y_var{this_var_pair})
                error('wrong variable name entered, no plot created');
            end
            ax(subplot_counter)= subplot(length(epochs),length(x_var),subplot_counter);
            hold(ax(subplot_counter),'on');
            
            % select cells
            if p.Results.common_cells_toggle % keep only cells common to selected epochs
                epoch_cells = intersect(remapping(epoch_idx(this_epoch),track_pair).new_ID,common_epoch_cells);
            elseif p.Results.session ~= 0 % look at single session
                session_idx = find(ismember(remapping(epoch_idx(this_epoch),track_pair).experiment,experiment_idx));
                epoch_cells= remapping(epoch_idx(this_epoch),track_pair).new_ID(session_idx);
            else % take all cells
                epoch_cells= remapping(epoch_idx(this_epoch),track_pair).new_ID;
                session_idx = 1:length(epoch_cells);
            end
            
            % Subsets
            if ~isempty(subset)
                    epoch_cells= intersect(epoch_cells,subset);
                    common_subset_cells= multintersect(remapping(epoch_idx(this_epoch),track_pair).new_ID(session_idx),epoch_cells,subset);    
            else % if no subset, just put cells common to epochs as default
                    common_subset_cells = intersect(remapping(epoch_idx(this_epoch),track_pair).new_ID(session_idx),epoch_cells);
            end
            
            % maybe not most straightforward way but makes sense
            epoch_cells_ind= session_idx(find(ismember(remapping(epoch_idx(this_epoch),track_pair).new_ID(session_idx),epoch_cells)));
            common_subset_cells_ind= session_idx(find(ismember(remapping(epoch_idx(this_epoch),track_pair).new_ID(session_idx),common_subset_cells)));
            
            if size(p.Results.cdata,1)==1
                cdata{epoch_idx(this_epoch)}= repmat(p.Results.cdata,[length(remapping(epoch_idx(this_epoch),track_pair).new_ID(session_idx)) 1]);
                cdata_epoch_idx = find(ismember(remapping(epoch_idx(this_epoch),track_pair).new_ID(session_idx),epoch_cells));
                cdata_common_subset_idx = find(ismember(remapping(epoch_idx(this_epoch),track_pair).new_ID(session_idx),common_subset_cells));
            else
                cdata= p.Results.cdata;
                cdata_epoch_idx = find(ismember(remapping(epoch_idx(this_epoch),track_pair).new_ID(session_idx),epoch_cells));
                cdata_common_subset_idx = find(ismember(remapping(epoch_idx(this_epoch),track_pair).new_ID(session_idx),common_subset_cells));
            end
            if ~isempty(p.Results.colormap)
                colormap(p.Results.colormap);
            end
            % plot all available cells as empty circles (selected as common to epochs or not based on inputs)
            scatter(ax(subplot_counter),...
                        remapping(epoch_idx(this_epoch),track_pair).(x_var{this_var_pair})(epoch_cells_ind), remapping(epoch_idx(this_epoch),track_pair).(y_var{this_var_pair})(epoch_cells_ind),...
                        20,cdata{epoch_idx(this_epoch)}(cdata_epoch_idx,:),...
                        'MarkerFaceAlpha',0,'MarkerEdgeAlpha',1);
           % plot selected cells as filled circles  
            scatter(ax(subplot_counter),...
                        remapping(epoch_idx(this_epoch),track_pair).(x_var{this_var_pair})(common_subset_cells_ind), remapping(epoch_idx(this_epoch),track_pair).(y_var{this_var_pair})(common_subset_cells_ind),...
                        20,cdata{epoch_idx(this_epoch)}(cdata_common_subset_idx,:),'filled',...
                        'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
            index_non_NaNs= intersect(find(~isnan(remapping(epoch_idx(this_epoch),track_pair).(x_var{this_var_pair}))),common_subset_cells_ind);
             if ~isempty(index_non_NaNs)
                    lm = fitlm(remapping(epoch_idx(this_epoch),track_pair).(x_var{this_var_pair})(index_non_NaNs),...
                        remapping(epoch_idx(this_epoch),track_pair).(y_var{this_var_pair})(index_non_NaNs),'linear');
                    %[pval(this_epoch,this_var_pair),F(this_epoch,this_var_pair),~] = coefTest(lm);
                    [pval(this_epoch,this_var_pair,track_pair),F(this_epoch,this_var_pair,track_pair),~] = coefTest(lm);
                    x=[min(remapping(epoch_idx(this_epoch),track_pair).(x_var{this_var_pair})(common_subset_cells_ind)) max(remapping(epoch_idx(this_epoch),track_pair).(x_var{this_var_pair})(common_subset_cells_ind))];
                    b=lm.Coefficients.Estimate';
                    plot(x,polyval(fliplr(b),x),'k--');
                    %slope_intcpt_lm(this_epoch,this_var_pair,:)= fliplr(b); 
                    slope_intcpt_lm(this_epoch,this_var_pair,:,track_pair)= fliplr(b);
%                     num_cells= lm.NumObservations
                    %title([epochs{this_epoch} ', pval= ' num2str(pval(this_epoch,this_var_pair),2)]);
                    title([epochs{this_epoch}]);
                    %text(ax(subplot_counter),.7,0.1,['p = ' num2str(pval(this_epoch,this_var_pair),2)],'Units','Normalized','FontSize',12,'FontName','Arial');
                    text(ax(subplot_counter),.7,0.1,['p = ' num2str(pval(this_epoch,this_var_pair,track_pair),2)],'Units','Normalized','FontSize',12,'FontName','Arial');

             else
                 title(epochs{this_epoch})
             end
             if isempty(p.Results.x_label)
                 xlabel_stripped= x_var{this_var_pair};
                 xlabel_stripped(xlabel_stripped == '_')= ' ';
                 xlabel(xlabel_stripped);
             else
                 if ischar(p.Results.x_label)
                     xlabel(p.Results.x_label);
                 else
                     xlabel(p.Results.x_label{this_var_pair});
                 end
             end
             if isempty(p.Results.y_label)
                 ylabel_stripped= y_var{this_var_pair};
                 ylabel_stripped(ylabel_stripped == '_')= ' ';
                 ylabel(ylabel_stripped);
             else
                 if ischar(p.Results.y_label)
                     ylabel(p.Results.y_label);
                 else
                     ylabel(p.Results.y_label{this_var_pair});
                 end
             end
             if p.Results.colorbar && subplot_counter==(length(epochs)*length(x_var))% last subplot
                 colorbar;
             end
             %xl= xlim; xticks([xl(1) 0 xl(2)]); 
             xtickformat('%i');
             %yl= ylim; yticks([yl(1) 0 yl(2)]); 
             ytickformat('%0.2f');

%              if strcmp(y_var{this_var_pair},'mean_max_FR_replay_diff')
% %                  yt=cellfun(@str2num,yticklabels).*sum(gausswin(0.25/(1/1000)));
% %                  yticklabels(yt); % this to account for bad nomalisation
% %                   ytickformat('%0.2f');
%              end
             
             %set(ax(subplot_counter),'TickDir','out','FontName','Arial','FontSize',12,'TickLength',[0.01 1]);
             %offsetAxes(ax(subplot_counter),5);
             subplot_counter= subplot_counter+1;
        end
    end
      run_format_settings(f1,'match_ax')
end

%run_format_settings(f1,'match_ax')
% if any(ismember(y_var,'mean_max_FR_replay_diff'))
%     allax = flipud(findall(f1,'type','axes'));
%     ib = find(ismember(y_var,'mean_max_FR_replay_diff'))*[1:length(epochs)];
%     lnk_ay = allax(ib);
%     %yl = arrayfun(@(x) allax(x).YLabel.String, 1:length(allax),'UniformOutput',0);
%     %[~,ib] = ismember(yl,'peak replay rate difference (Hz)'); % this should be the final Y label we settle on
%     %lnk_ay = allax(logical(ib));
%     for j = 1 : length(lnk_ay)
%         yt = cellfun(@str2num,yticklabels(lnk_ay(j))).*sum(gausswin(0.25/(1/1000)));
%         lnk_ay(j).YTickLabels = round(yt,2);
% %         yticklabels(yt); % this to account for bad nomalisation
%     end
% end

end