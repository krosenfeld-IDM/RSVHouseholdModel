function [Q] = prediction_interval_for_mixed_poisson(spread,lower_Q,upper_Q)
%Given an empirical spread of predictions for mean value of a Poisson I
%use the empirical CDF to calculate the lower and upper
%quantiles

f_lower = @(x) empirical_cdf_possions(spread,x) - lower_Q;
lower = fzero(f_lower,10000);

f_upper = @(x) empirical_cdf_possions(spread,x) - upper_Q;
upper = fzero(f_upper,10000);

Q = [lower,upper];

end

