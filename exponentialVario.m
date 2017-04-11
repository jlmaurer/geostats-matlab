function [var, cov] = exponentialVario(param, h, include_nugget)
% This funparamtion returns the variogram and covariance for the exponential
% variogram. 

if nargin<3, include_nugget = 0; end

if include_nugget==1
    var = param(1)*(1-exp(-h/param(2))) + param(3)*(h~=0);
    cov = param(1)*exp(-h/param(2)) + param(3)*(h==0); 
else
    var = param(1)*(1-exp(-h/param(2)));
    cov = param(1)*exp(-h/param(2)); 
end

end
