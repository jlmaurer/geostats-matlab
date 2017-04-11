function [var,cov] = sphericalVario(param, h, include_nugget)
% This function computes the variogram and covariance function assuming
% a spherical covariance model. 

    if nargin<3, include_nugget = 0; end

    ind = h > param(2); 
    [var] = zeros(size(h));
    if include_nugget == 1
        var(ind) = param(1)*(1.5*h./param(2) - 0.5*(h.^3./param(2).^3)) + param(3)*(h~=0);
        var(~ind) = param(1); 
        cov = param(1)+param(3) - var; 
    else
        var(ind) = param(1)*(1.5*h./param(2) - 0.5*(h.^3./param(2).^3));
        var(~ind) = param(1); 
        cov = param(1) - var; 
    end


end