function [hEff,rawVario, trend] = rawvario(X,data, MT,trend)
% raw semivariances of the data located at X.
%
% Inputs: 
%   X       - an MxD matrix, where M is the number of observations
%               and D is the dimension of space
%   data    - Mx1 vector of values at each observation location
%   MT      - an anisotropy matrix to be applied to X
%   trend   - a structure containing fields: 
%               flag:   0/1         indicates whether there is a trend
%               xy:     same as X   Locations to estimate trend if flag =1 
%               At this time only linear trends are supported
%
% Outputs: 
%   hEff    - effective lag distances between all observations
%   varVario- raw semivariances between all data
%   trend   - same structure as input but with the additional fields:
%               Trend_params:   parameter vector for the trend estimate
%               Dpred:          predicted data (due to just the trend) at each X
%
%REVISION HISTORY:
%   tae, 3/10/99 created. --> original code was obtained in the Geostatistics class
%       taught by Anna Michalak at Stanford University
%   JLM, 4/11/2017 modified  and name changed from varioraw_

% set defaults
if nargin<4,trend = struct('flag',0); end
if nargin<3, MT=[]; end

% subtract trend if needed
if trend.flag ==1
    G= [ones(size(trend.xy,1), 1), trend.xy(:,1), trend.xy(:,2)];
    beta = G\data; 
    dpred = G*beta; 
    dataDetrend= data - dpred; 
    trend.Trend_Params=beta;
    trend.Dpred = dpred;
else
    dataDetrend=data;
end

% build raw variogram
rawVario = 0.5*distance_(dataDetrend,dataDetrend).^2;

% if anisotropy matrix provided, transform to isotropic coordinates
if ~isempty(MT)
    XEff=(MT*X')';
else
    XEff = X; 
end
hEff = distance_(XEff,XEff);

end
