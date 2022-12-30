

function form = formatting_settings

%%% MAIN FIGURE
% Axis settings
form.TickDir = 'out';
form.TickLength = [.005 1];
form.LineWidth = 1.5;
form.tickformat_dec = '%0.2f';
form.tickformat = '%i';
form.NumYTicks = 3;
form.NumXTicks = 5; %5

% Font settings
form.FontSize = 12;
form.FontName = 'Arial';

%%% OTHER 
% add the type of figure as a field

form.ColorbarTicks = 3;

% Inset properties
form.inset.FontSize = 11;
form.inset.NumYTicks = 3;
form.inset.NumXTicks = 5;
form.inset.TickDir = 'out';
form.inset.TickLength = [.005 1];
form.inset.LineWidth = 1.5;

end
