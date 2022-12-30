function [value_at_p,p_idx]= get_percentiles_pdf(pd,centres,p)

cum_pd= cumsum(pd);
p_score= sum(pd)*p/100;
[~,p_idx]= min(abs(cum_pd - p_score));
value_at_p= pd(p_idx);
p_idx= centres(p_idx);

end