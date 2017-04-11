function [var,cov] = maternVario(param,nu, h, include_nugget)
% This function returns the variogram and covariance at h
% for the Matern covariance model with parameters param. If
% include_nugget is 1, a nugget is included. 
%
% param should be in the form [sill, range, nu parameter, nugget]

if nargin < 4, include_nugget = 0; end

factor1 = - (nu-1)*log(2) - log(gamma(nu)); 
factor2 = nu*(log(h) - log(param(2))); 
factor3 = log(besselk(nu, (h./param(2)))); 

term = factor1 + factor2 + factor3; 

V1 = param(1)*(1 - exp(term));
V1(isnan(V1)) = 0; 
C1 = param(1) - V1; 

if include_nugget==1
    var = V1 + param(3)*(h~=0); 
    cov = C1 + param(3)*(h==0);
else
    var = V1; 
    cov = C1; 
end


end
