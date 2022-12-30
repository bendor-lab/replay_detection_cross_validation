% LINE FITTING
% Find best-fit-line using radon transform.
    % The output will be a 'R' matrix where each column is the radon transform per degree. 
    % and 'XP', which are the radian coordinates along the x axis.
% Loads 'decoded_replay_events'
% Adds line fit score for each replay event to a new column in replay_track structure

function [linear_score, best_line_pixels] = line_fitting(decoded_event)
         linear_score=NaN;
        %for k=1:size(decoded_event,2) %number of time bins
          % decoded_event_smooth(:,k)=3*smooth(decoded_event(:,k),3); %each bin is the sum of the current bin plus the previous and next one
        %end
        gsimage = mat2gray((decoded_event)); %converts matrix to grayscale image with values from 0 to 1
        theta = 1:1:179; % changed from 0.1 steps so that it can run faster . avoid horizontal lines (theta = 0 or 180)
        [R,xp] = radon(gsimage,theta); 
        
        % Finds maximum of radon transform- the line of best fit is orthogonal to this angle.
        % This line is centred in the middle of the plot- need to recalculate y intercept
        
        [M, I] = max(R(:)); %creates a column from the R matrix and finds maximum
        [X, T] = ind2sub(size(R), I);%finds indices of max value back in the R matrix
        offset = xp(X);  %offset from midpoint along slice (radian in max point)
        slope = cotd(theta(T)); %in pixels
        y0 = (floor(size(decoded_event,1) + 1) /2) - offset * sind(theta(T)) - (floor(size(decoded_event,2) + 1) /2 + offset * cosd(theta(T))) * slope;
       
        best_line_pixels = [slope, y0];  %in pixels, not position/time 
        
        % output = M/size(decoded_event_smooth,2);
        
          %find the closest bin and sum probability
       y1=floor(polyval([slope y0],1:size(decoded_event,2)));
       y2=ceil(polyval([slope y0],1:size(decoded_event,2)));
       y=[y1 y2]; %up to 1 neighboring bins (within 10-20cm of location)
       x=[1:size(decoded_event,2) 1:size(decoded_event,2)];

      index=find(y>size(decoded_event,1)  | y<1 | isnan(y));
      y(index)=[];
      x(index)=[];
      decoding_mask=zeros(size(decoded_event));
      for c=1:length(y)
          decoding_mask(y(c),x(c))=1;
      end
      linear_score=sum(sum(decoded_event.*decoding_mask))/size(decoded_event,2);
        
        
end