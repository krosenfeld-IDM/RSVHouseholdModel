function [P] = empirical_cdf_possions(spread,x)
%Given an empirical spread of predictions for mean value of a Poisson I
%construct the empirical CDF
P = mean(poisscdf(x*ones(size(spread)),spread));

end

