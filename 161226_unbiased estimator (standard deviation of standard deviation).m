clear all

%https://en.wikipedia.org/wiki/Unbiased_estimation_of_standard_deviation

 n=4;    %the number used for 4-point calc
 Sds_coef=(1-2/(n-1)*(gamma(n/2)/gamma((n-1)/2))^2)^0.5

 Sds_unbaised_coef=gamma((n-1)/2)/gamma(n/2)*((n-1)/2-(gamma(n/2)/gamma((n-1)/2))^2)^0.5