function [Dest, Dsig, lambda] = kriging_cpuk(SIG, sig0, d, sig2,xy, XY, dim)
% This function returns weights to use for continous-part Kriging
% w/variable mean, using the covariance matrices SIG (covariance of the
% observation locations) and sig0 (covariance of the observation locations
% with the measurement locations). All of the estimation locations XY (a 
% Nestx2 matrix) are computed simultaneously. At this time, only linear 
% trends are supported by the code. 
%
% Inputs:
%   SIG     - Matrix of correlations between observation locations (if
%               there are M observed sites, this is MxM)
%   sig0    - matrix of correlations between observed locations and
%               estimation locations. If there are N estimation locations, 
%               this is an MxN matrix. 
%   d       - Observed values (Mx1)
%   sig2    - Point-wise variance, or covariance function at distance 0
%   xy      - observation locations 
%   XY      - estimation locations
%   dim     - if dim is specified, it is the xy-component to estimate the
%               linear trend (can be 1 or 2). If not specified, a trend in 
%               both component directions is estimated. 
%
% Outputs: 
%   Dest    - Estimated values at the N estimation locations (Nx1)
%   Dsig    - Estimated uncertainty (1-sigma) of Dest
%   lambda  - matrix of weights used to compute the estimated values
%
% Author: Jeremy Maurer, April 7, 2017
% License: MIT

    % specify parameters
    Nobs = size(SIG,2); 
    Nest = size(sig0, 2); 
    
    if nargin < 7
        dim=2;
    end
    
    % check which directions to estimate a linear trend
    if dim==2
        F = [ones(Nobs,1), xy(:,1), xy(:,2)]; 
        f0 = [ones(1,Nest); XY(:,1)'; XY(:,2)']; 
    elseif dim == 1
        F = [ones(Nobs,1), xy(:,1)]; 
        f0 = [ones(1,Nest); XY(:,1)']; 
    else
        F =  [ones(Nobs,1), xy(:,2)]; 
        f0 = [ones(1,Nest); XY(:,2)']; 
    end
    
    % build matrices
    A= [SIG, F; F', zeros(size(F,2))]; 
    B = [sig0; f0]; 

    % solve for weights
    lambdanu = A\B; 
    lambda = lambdanu(1:size(sig0,1),:); 
    nu = -lambdanu(size(sig0,1)+1:end,:); 
    fnu = sum(f0.*nu, 1); 
    
    % compute estimated data and uncertainty
    Dest = lambda'*d; 
    Dsig = sqrt(sig2 - diag(lambda'*sig0) + fnu(:)); 
end