function [var, cov] = linearVario(param,h, include_nugget)
% this function computes a linear semivariance function

    if nargin < 3, include_nugget = 0; end
    
    [var] = deal(zeros(size(h))); 
    if include_nugget ==1 
        ind = h<param(2); 
        hind = h(ind); 
        var(ind) = (param(1)./param(2))*hind + param(3)*(hind~=0);
        var(~ind) = param(1)+ param(3);
        cov = param(3) + param(1) - var;
    else
        ind = h<param(2); 
        var(ind) = (param(1)./param(2))*h(ind);
        var(~ind) = param(1); 
        cov = param(1) - var;
    end

end