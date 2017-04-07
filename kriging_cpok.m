function [Dest, Dsig, lambda] = kriging_cpok(SIG, sig0, d, sig2)
% This function returns weights to use for continous-part ordinary Kriging, 
% using the covariance matrices SIG (covariance of the
% observation locations) and sig0 (covariance of the observation locations
% with the measurement locations).All of the estimation locations are
% computed simultaneously. 
%
% Inputs:
%   SIG     - Matrix of correlations between observation locations (if
%               there are M observed sites, this is MxM)
%   sig0    - matrix of correlations between observed locations and
%               estimation locations. If there are N estimation locations, 
%               this is an MxN matrix. 
%   d       - Observed values (Mx1)
%   sig2    - Point-wise variance, or covariance function at distance 0
%
% Outputs: 
%   Dest    - Estimated values at the N estimation locations (Nx1)
%   Dsig    - Estimated uncertainty (1-sigma) of Dest
%   lambda  - matrix of weights used to compute the estimated values
%
% Author: Jeremy Maurer, April 7, 2017
% License: MIT

    Nobs = size(SIG,2); 
    Nest = size(sig0, 2); 

    A= [SIG, ones(Nobs,1); ones(1,Nobs), 0]; 
    B = [sig0; ones(1,Nest)]; 

    lambdanu = A\B; 
    lambda = lambdanu(1:end-1,:); 
    nu = -lambdanu(end,:); 
    
    % compute predicted velocity and uncertainty
    Dest = lambda'*d; 
    Dsig = sqrt(sig2 - diag(lambda'*sig0) + nu(:)); 
end