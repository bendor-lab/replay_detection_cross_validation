

function run_format_settings(figure_handle,varargin)
%varargin - settings that are not for the main figure
% match_ax, match_Xax, match_Yax: match limits of multiple axes. Can choose both, or only X or Y
% ...any other type of figure that doesn't run with the main settings

%%%% Load settings
form = formatting_settings;
fieldsform = fields(form);
if isempty(figure_handle)
    f = gcf;
else
    f = figure_handle;
end
allax = findall(f,'type','axes');

% Find if there are insets
inset_ax= [];
for k = 1 : length(allax) %find axis indx of inset
    idx = find(1: length(allax) ~= k);
    inset_ax = [inset_ax idx(cell2mat(arrayfun(@(x)  allax(k).Position(1) < allax(x).Position(1) &  allax(k).Position(2) < allax(x).Position(2) &...
        (allax(k).Position(1)+allax(k).Position(3)) >= (allax(x).Position(1)+allax(x).Position(3)),idx,'UniformOutput',0)))];
    % &  (allax(k).Position(2)+allax(k).Position(4)) <=  (allax(x).Position(2)+allax(x).Position(4))
end
if ~isempty(inset_ax)
    main_ax = find(~ismember(1: length(allax),unique(inset_ax))); %exclude axis indx of inset
else
    main_ax = 1: length(allax);
end

%%%%% General figure settings
arrayfun(@(x) set(allax(x),'box','off'),1:length(allax),'UniformOutput',0);
colB_idx = find(cell2mat(arrayfun(@(x) ~isempty(allax(x).Colorbar),1:length(allax),'UniformOutput',0)));
if ~isempty(colB_idx) % change colorbar settings, if exists
    for cb =  1 : length(colB_idx)
        allax(colB_idx(cb)).Colorbar.Ticks = linspace(allax(colB_idx(cb)).Colorbar.Ticks(1),allax(colB_idx(cb)).Colorbar.Ticks(end),form.ColorbarTicks);
        allax(colB_idx(cb)).Colorbar.TickLabels = num2cell(allax(colB_idx(cb)).Colorbar.Ticks);
        allax(colB_idx(cb)).Colorbar.FontSize = form.FontSize;
        allax(colB_idx(cb)).Colorbar.FontName = form.FontName;
    end
end


if ~isempty(varargin)
    for ax = 1 : length(allax)
        % Find indices of axes that share X or Y labels
        xl = arrayfun(@(x) allax(x).XLabel.String, main_ax,'UniformOutput',0);
        [~,ia] = ismember(xl,allax(ax).XLabel.String);
        lnk_ax = main_ax(logical(ia));
        yl = arrayfun(@(x) allax(x).YLabel.String, main_ax,'UniformOutput',0);
        [~,ib] = ismember(yl,allax(ax).YLabel.String);
        lnk_ay = main_ax(logical(ib));
        
        if length(allax) > 1 & strcmp(varargin,'match_ax') |  strcmp(varargin,'match_Xax') % match X axes limits
            lims = [allax(lnk_ax).XLim];
            [~,mxid] = max(lims);
            [~,mnid] = min(lims);
            [allax(lnk_ax).XLim] = deal([lims(mnid) lims(mxid)]);
        end
        if length(allax) > 1 & strcmp(varargin,'match_ax') |  strcmp(varargin,'match_Yax') % match Y axes limits
            multiple_yaxis_idx = cell2mat(arrayfun(@(x) contains(allax(x).YAxisLocation,'right'),lnk_ay,'UniformOutput',0));
            switch  any(multiple_yaxis_idx)
                case 0 %no multiple axes
                    lims = [allax(lnk_ay).YLim];
                    [~,mxid] = max(lims);
                    [~,mnid] = min(lims);
                    [allax(lnk_ay).YLim] = deal([lims(mnid) lims(mxid)]);
                case 1 %multiple axes
                    multiple_yaxis_idx = find(cell2mat(arrayfun(@(x) contains(allax(x).YAxisLocation,'right'),lnk_ay,'UniformOutput',0)));
                    single_yaxis_idx = find(~ismember(1: length(lnk_ay),unique(multiple_yaxis_idx)));
                    if ~isempty(single_yaxis_idx)
                        lims = [allax(single_yaxis_idx).YLim];
                        [~,mxid] = max(lims);
                        [~,mnid] = min(lims);
                        [allax(lnk_ay).YLim] = deal([lims(mnid) lims(mxid)]);
                    end
                    lims = [allax(multiple_yaxis_idx).YLim];
                    [~,mxid] = max(lims);
                    [~,mnid] = min(lims);
                    [allax(multiple_yaxis_idx).YLim] = deal([lims(mnid) lims(mxid)]);
                    arrayfun(@(x) yyaxis(allax(x),'left'),multiple_yaxis_idx,'UniformOutput',0) %change to left yaxis
                    lims = [allax(multiple_yaxis_idx).YLim];
                    [~,mxid] = max(lims);
                    [~,mnid] = min(lims);
                    [allax(multiple_yaxis_idx).YLim] = deal([lims(mnid) lims(mxid)]);
            end
        end
    end
end


%%%% Apply axis settings

%Find settings that refer to axis properties
prop_idx = find(cell2mat(arrayfun(@(x) isprop(allax(1),fieldsform{x}), 1:length(fieldsform),'UniformOutput',0)));

for j = 1 : length(allax)
    ax = allax(j);
    if any(ax == allax(main_ax))
        for k = 1 : length(prop_idx)
            set(ax,fieldsform{prop_idx(k)},form.(sprintf('%s',fieldsform{prop_idx(k)})))
        end
        num_yticks = form.NumYTicks;
        num_xticks = form.NumXTicks;
        
    elseif any(ax == allax(inset_ax)) %if inset
        
        subf = fields(form.inset);
        subf_prop = find(cell2mat(arrayfun(@(x) isprop(ax,subf{x}), 1:length(subf),'UniformOutput',0)));
        subf = subf(subf_prop);
        for j = 1 : length(subf)
            set(ax,subf{j},form.inset.(sprintf('%s',subf{j})))
            xl = get(ax,'xlabel');
            yl = get(ax,'ylabel');
        end
        num_yticks = form.inset.NumYTicks;
        num_xticks = form.inset.NumXTicks;
        Yax_id = ax;
        Xax_id = ax;
    end
    
    %%% Adjust axes ticks
    if length(ax.YTick) > num_yticks
        new_Yticks = linspace(min(ylim(ax)),max(ylim(ax)),num_yticks);
        if new_Yticks(2) < 0 & any(new_Yticks~=0) %if Y axis starts negative, make sure a tick is zero
            [~,mnid] = min(abs(new_Yticks));
            new_spacing = (0 - new_Yticks(1))/(mnid-1);
            new_ticks_sp = new_Yticks(1):new_spacing:new_Yticks(end);
            ax.YTick = new_ticks_sp;
        else
            ax.YTick = new_Yticks;
        end
    end
    
    if length(ax.XTick) > num_xticks
        new_Xticks = linspace(min(xlim(ax)),max(xlim(ax)),num_xticks);
        if new_Xticks(2) < 0 | any(new_Xticks~=0) %if X axis starts negative, make sure a tick is zero
            [~,mnid] = min(abs(new_Xticks));
            new_spacing = (0 - new_Xticks(1))/(mnid-1);
            new_ticks_sp = new_Xticks(1):new_spacing:max(xlim(ax));%new_Xticks(end);
            ax.XTick = new_ticks_sp;
        else
            ax.XTick = new_Xticks;
        end
    end
    
end
end

