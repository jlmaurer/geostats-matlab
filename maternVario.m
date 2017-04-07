function [var,cov] = maternVario(param, h, include_nugget)
% This function returns the variogram and covariance at h
% for the Matern covariance model with parameters param. If
% include_nugget is 1, a nugget is included. 
%
% param should be in the form [sill, range, nu parameter, nugget]

if nargin < 3, include_nugget = 0; end

factor1 = 1/((2^(param(3)-1))*gamma(param(3))); 
factor2 = (h./param(2)).^param(3); 
factor3 = besselk(param(3), (h./param(2))); 

V1 = param(1)*(1 - factor1*factor2.*factor3);
C1 = param(1) - V1; 

if include_nugget==1
    var = V1 + param(4)*(h~=0); 
    cov = C1 + param(4)*(h==1);
else
    var = V1; 
    cov = C1; 
end
end
