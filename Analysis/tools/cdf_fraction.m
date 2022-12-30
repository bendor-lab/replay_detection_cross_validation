function y_distr= cdf_fraction(pd,x_bins)
% values
% bins
for this_bin=1:length(x_bins)
    idx= find(pd < x_bins(this_bin));
    y_distr(this_bin)= length(idx)/length(pd);
end
    
end