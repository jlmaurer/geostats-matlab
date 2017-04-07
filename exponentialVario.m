function [var, cov] = exponentialVario(param, h)
% This funparamtion returns the variogram and covariance for the exponential
% variogram. 

var = param(1)*(1-exp(-h/param(2))) + param(3)*(h~=0);
cov = param(1)*exp(-h/param(2)) + param(3)*(h==0); 
end
