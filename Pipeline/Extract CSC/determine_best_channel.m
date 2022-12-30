function best_channels = determine_best_channel(varargin)
% determines best channel for theta, gamma, etc. based on the peak response of the PSD (relative to nearby channels)
%INPUTS:
%   - 'all','hpc', 'ofc' or 'vc'

load PSD_data
parameters= list_of_parameters;

if ~isempty(varargin)
    switch varargin{1}
        case 'all'
            MODE = [1 1];
        case 'hpc'
            MODE = [1 0];
        case 'ofc'
            MODE = [0 1];
        case 'vc'
            MODE = [0 1];
    end
else
    MODE= [1 1];
end

% Looks into first CSC frequencies and finds periods with theta, ripple and delta oscillations
if MODE(1)==1 %analyze theta, ripple and delta
    if isfield(all_PSD,'PSD_F_awake')
        theta_index=find(all_PSD(1).PSD_F_awake>=parameters.theta_filter(1) & all_PSD(1).PSD_F_awake<=parameters.theta_filter(2));
        ripple_index=find(all_PSD(1).PSD_F_sleep>=parameters.ripple_filter(1) & all_PSD(1).PSD_F_sleep<=parameters.ripple_filter(2));
        delta_index=find(all_PSD(1).PSD_F_sleep>=parameters.delta_filter(1) & all_PSD(1).PSD_F_sleep<=parameters.delta_filter(2));
    else
        theta_index=find(all_PSD(1).PSD_F >=parameters.theta_filter(1) & all_PSD(1).PSD_F <=parameters.theta_filter(2));
        ripple_index=find(all_PSD(1).PSD_F >=parameters.ripple_filter(1) & all_PSD(1).PSD_F <=parameters.ripple_filter(2));
        delta_index=find(all_PSD(1).PSD_F >=parameters.delta_filter(1) & all_PSD(1).PSD_F <=parameters.delta_filter(2));
    end
end
if MODE(2)==1 %analyze spindle and high gamma
    if isfield(all_PSD,'PSD_F_awake')
        spindle_index=find(all_PSD(1).PSD_F_sleep>=parameters.spindle_filter(1) & all_PSD(1).PSD_F_sleep<=parameters.spindle_filter(2));
        high_gamma_index=find(all_PSD(1).PSD_F_awake>=parameters.high_gamma_filter(1) & all_PSD(1).PSD_F_awake<=parameters.high_gamma_filter(1));
    else
        spindle_index=find(all_PSD(1).PSD_F>=parameters.spindle_filter(1) & all_PSD(1).PSD_F<=parameters.spindle_filter(2));
        high_gamma_index=find(all_PSD(1).PSD_F>=parameters.high_gamma_filter(1) & all_PSD(1).PSD_F<=parameters.high_gamma_filter(1));
    end 
end

% For each type of oscillation, calculate peak of power. 

% this is achieved by comparing the median overall power within the frequency range
% with the mean of only the lowest and highest channel of the frequency range. 
% If there is an increase of power, relative to the 1/f decrease in power over the frequency range, the
% overall median will be higher than the mean of the first and last values

for thisCSC=1:length(all_PSD)
    if MODE(1)==1  %analyze theta, ripple and delta
        if isfield(all_PSD,'PSD_F_awake')
            % choose theta channel based on awake PSD, but delta and ripple based on sleep PSD 
            y.theta(thisCSC,:) = median(all_PSD(thisCSC).PSD_awake(:,theta_index),2)-0.5*(all_PSD(thisCSC).PSD_awake(:,theta_index(1))+all_PSD(thisCSC).PSD_awake(:,theta_index(end)));
            y.ripple(thisCSC,:)= median(all_PSD(thisCSC).PSD_sleep(:,ripple_index),2)-0.5*(all_PSD(thisCSC).PSD_sleep(:,ripple_index(1))+all_PSD(thisCSC).PSD_sleep(:,ripple_index(end)));
            y.delta(thisCSC,:) = median(all_PSD(thisCSC).PSD_sleep(:,delta_index),2)-0.5*(all_PSD(thisCSC).PSD_sleep(:,delta_index(1))+all_PSD(thisCSC).PSD_sleep(:,delta_index(end)));
        else
            y.theta(thisCSC,:) = median(all_PSD(thisCSC).PSD(:,theta_index),2)-0.5*(all_PSD(thisCSC).PSD(:,theta_index(1))+all_PSD(thisCSC).PSD(:,theta_index(end)));
            y.ripple(thisCSC,:)= median(all_PSD(thisCSC).PSD(:,ripple_index),2)-0.5*(all_PSD(thisCSC).PSD(:,ripple_index(1))+all_PSD(thisCSC).PSD(:,ripple_index(end)));
            y.delta(thisCSC,:) = median(all_PSD(thisCSC).PSD(:,delta_index),2)-0.5*(all_PSD(thisCSC).PSD(:,delta_index(1))+all_PSD(thisCSC).PSD(:,delta_index(end)));
        end
        
        y.theta_power(thisCSC,:) = mean(maxk(y.theta(thisCSC,:),4));
        y.ripple_power(thisCSC,:)= mean(maxk(y.ripple(thisCSC,:),4));
        y.delta_power(thisCSC,:) = mean(maxk(y.delta(thisCSC,:),4));
    end
    if MODE(2)==1 %analyze spindle and high gamma
        if isfield(all_PSD,'PSD_F_awake')
            y.spindle(thisCSC,:)=median(all_PSD(thisCSC).PSD_sleep(:,spindle_index),2)-0.5*(all_PSD(thisCSC).PSD_sleep(:,spindle_index(1))+all_PSD(thisCSC).PSD_sleep(:,spindle_index(end)));
            y.high_gamma(thisCSC,:)=median(all_PSD(thisCSC).PSD_awake(:,high_gamma_index),2)-0.5*(all_PSD(thisCSC).PSD_awake(:,high_gamma_index(1))+all_PSD(thisCSC).PSD_awake(:,high_gamma_index(end)));
        else
            y.spindle(thisCSC,:)=median(all_PSD(thisCSC).PSD(:,spindle_index),2)-0.5*(all_PSD(thisCSC).PSD(:,spindle_index(1))+all_PSD(thisCSC).PSD(:,spindle_index(end)));
            y.high_gamma(thisCSC,:)=median(all_PSD(thisCSC).PSD(:,high_gamma_index),2)-0.5*(all_PSD(thisCSC).PSD(:,high_gamma_index(1))+all_PSD(thisCSC).PSD(:,high_gamma_index(end)));
        end
        y.spindle_power(thisCSC,:)=mean(maxk(y.spindle(thisCSC,:),4));
        y.high_gamma_power(thisCSC,:)=mean(maxk(y.high_gamma(thisCSC,:),4));
    end
end


if MODE(1)==1
    % Find the CSC with the highest power for each oscillation
    [~,y.bestCSC_theta]  = max(y.theta_power);
    [~,y.bestCSC_ripple] = max(y.ripple_power);
    [~,y.bestCSC_delta]  = max(y.delta_power);
    
    % Find CSC with highest theta and lowest ripple power (used for replay event detection)
    [~,min_index] = sort(y.ripple_power,'ascend'); %sort from min ripple to max ripple power
    [~,max_index] = sort(y.theta_power,'descend'); %sort from max theta to min theta power
    
    all_positions = [];
    for CSC=1:length(min_index)
        CSC_position = (find(min_index == CSC)) + (find(max_index == CSC));
        all_positions = [all_positions, CSC_position];
    end
    [~,y.bestCSC_ThetaRippleDiff_RankMethod] = min(all_positions);
    
    %calculate channel with biggest difference between theta and ripple normalized power.
     y.normalized_theta_power=y.theta_power/max(y.theta_power);  %power normalized by max power across channels
     y.normalized_ripple_power=y.ripple_power/max(y.ripple_power);  %power normalized by max power across channels
      y.normalized_ripple_power(find(y.normalized_ripple_power<0))=0; %if normalized ripple power is less than zero, treat it as 0
     [~,y.bestCSC_ThetaRippleDiff_NormalizedMethod]  = max(y.normalized_theta_power-y.normalized_ripple_power);  %find channel with maximum power difference between theta and ripple
     
    %convert index to channel
    y.bestCSC_theta=all_PSD(y.bestCSC_theta).CSCchannel;
    y.bestCSC_ripple=all_PSD(y.bestCSC_ripple).CSCchannel;
    y.bestCSC_delta=all_PSD(y.bestCSC_delta).CSCchannel;
    y.bestCSC_ThetaRippleDiff_RankMethod=all_PSD(y.bestCSC_ThetaRippleDiff_RankMethod).CSCchannel;
     y.bestCSC_ThetaRippleDiff_NormalizedMethod=all_PSD(y.bestCSC_ThetaRippleDiff_NormalizedMethod).CSCchannel;
     
elseif MODE(2)==1
    % Find the CSC with the highest power for each oscillation
    [~,y.bestCSC_spindle] = max(y.spindle_power);
    [~,y.bestCSC_high_gamma] = max(y.high_gamma_power);
    
    %convert index to channel
    y.bestCSC_spindle=all_PSD(y.bestCSC_spindle).CSCchannel;
    y.bestCSC_high_gamma=all_PSD(y.bestCSC_high_gamma).CSCchannel;
    
end

% save as
best_channels = y;
save('best_CSC','best_channels');

end