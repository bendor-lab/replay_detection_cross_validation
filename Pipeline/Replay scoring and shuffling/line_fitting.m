function [linear_score,out]=line_fitting(decoded_event)

FIGURES_ON=0;
[position_bins,time_bins]=size(decoded_event);
time_bins2 = min([time_bins position_bins]);

max_value=max([time_bins2 position_bins]);
slope=(2/max_value):(2/max_value):.999;

slope=[slope 1 fliplr(1./slope)];


if FIGURES_ON
    figure
end
for s=1:length(slope)
    kernel=zeros(position_bins,time_bins2);
    x=1:time_bins2;
    y=round(slope(s)*(x-1))+1;
    %I rounded the line to the nearest pixel
    k=[];
    for i=1:time_bins2
        j=y(i);
        if j<=position_bins
            kernel(j,i)=1;   %NEED TO PUT J BEFORE I???
            if k<j
                kernel(k+1:j,i)=1;
            end
            k=j;
        elseif k~=position_bins
            j=position_bins;
            kernel(j,i)=1;
            if k<j
                kernel(k+1:j,i)=1;
            end
            k=j;
        end
    end
    
    if FIGURES_ON
        imagesc(kernel); axis xy; pause(0.2);
    end
    %by changing line_width you could include the pixels above and below
    filter_result_P=filter2(kernel,decoded_event,'full');
    [score_P(s),index_P(s)]=max(filter_result_P(:));  %positive slope
    [I,J]=ind2sub(size(filter_result_P),index_P(s));
    y0_P(s)=-position_bins + I + (time_bins2 - J)*slope(s);
    
    filter_result_N = filter2(flipud(kernel), decoded_event,'full');
    [score_N(s), index_N(s)]=max(filter_result_N(:));  %negative slope
    [I,J]=ind2sub(size(filter_result_N),index_N(s));
    y0_N(s)= I -(time_bins2 - J)*slope(s);
end



if max(score_P)>=max(score_N)
    [linear_score,index]=max(score_P);
    out.slope = slope(index);
    out.y0 = y0_P(index);
    if FIGURES_ON
        figure
        imagesc(decoded_event);
        hold on
        plot([0.5,time_bins + 0.5],[out.y0 + 0.5,out.y0 + out.slope*time_bins + 0.5]); %everything is off by half a pixel
    end
else
    [linear_score,index]=max(score_N);
    out.slope = -slope(index);
    out.y0 = y0_N(index);
    if FIGURES_ON
        figure
        imagesc(decoded_event);
        hold on
        plot([0.5,time_bins + 0.5],[out.y0 + 0.5,out.y0 + out.slope*time_bins + 0.5]);
    end
end


end
