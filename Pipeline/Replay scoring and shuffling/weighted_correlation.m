% WEIGHTED CORRELATION
% Linear correlation between time and position weighted by the associated posterior probability
% INPUTS: decoded_position (e.g. decoded replay event); T as vector of time bins; P as vector of position bins

function output = weighted_correlation(decoded_event)

pixel_size = size(decoded_event);
P = 1:pixel_size(1); % position bins
T = 1:pixel_size(2); % time bins
decoded_event_sum=sum(sum(decoded_event));
% weighted mean of decoded time
for tbin = 1:pixel_size(2)
    probT(:,tbin) = decoded_event(:,tbin).*T(tbin);
end
weighted_mean_time =  sum(sum(probT))./decoded_event_sum;

% weighted mean of decoded position
for xbin = 1:pixel_size(1)
    probP(xbin,:) = decoded_event(xbin,:).*P(xbin);
end
 weighted_mean_position = sum(sum(probP))./decoded_event_sum;
 
%weighted covariance between time and decoded position
for tbin = 1:pixel_size(2)
    for xbin = 1:pixel_size(1)
        bin = decoded_event(xbin,tbin);
        probTP(xbin,tbin) = bin*(T(tbin)-weighted_mean_time).*(P(xbin)-weighted_mean_position);
    end
end
 cov_TP = sum(sum(probTP))./decoded_event_sum;

 % weighted covariance between decoded position and decoded position
for tbin = 1:pixel_size(2)
    for xbin = 1:pixel_size(1)
        bin = decoded_event(xbin,tbin);
        probPP(xbin,tbin) = bin*(P(xbin)-weighted_mean_position)^2;
    end
end
   cov_PP = sum(sum(probPP))./decoded_event_sum;
 
%weighted covariance between time and time
for tbin = 1:pixel_size(2)
    for xbin = 1:pixel_size(1)
        bin = decoded_event(xbin,tbin);
        probTT(xbin,tbin) = bin*(T(tbin)-weighted_mean_time)^2;
    end
end
 cov_TT = sum(sum(probTT))./decoded_event_sum;
 
% weighted correlation between time and decoded position
output = abs(cov_TP / sqrt(cov_TT*cov_PP));
end


