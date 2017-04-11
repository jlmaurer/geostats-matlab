function [var, K] = powerVario(param,h)
% this function computes a linear semivariance function

    var = param(1)*(h.^param(2));
    K = -var;

end