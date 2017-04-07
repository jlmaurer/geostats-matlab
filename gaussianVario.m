function [var, cov] = gaussianVario(param, h, thresh,include_nugget)
% This function returns the gaussian variogram and covariance for the
% specified parameters.

if nargin < 4
    if nargin < 3
        thresh = 0; 
    end
    include_nugget = 1; 
end

if include_nugget == 1
    var = param(1)*(1-exp(-(h.^2)./(param(2)^2))) + param(3)*(h~=0);
    cov = param(1)*exp(-(h.^2)/(param(2)^2)) + param(3)*(h==0);
else
    var = param(1)*(1-exp(-(h.^2)/(param(2)^2)));
    cov = param(1)*exp(-(h.^2)/(param(2)^2));
end

cov(cov < thresh) = 0; 

end