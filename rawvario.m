function [hEff,rawVario, trend] = rawvario(yLoc,yVal, MT,trend)

%REVISION HISTORY:
%   tae, 3/10/99 created. 
%   JLM, 4/7/2017 modified  

% set defaults
if nargin<4,trend = struct('flag',0); end
if nargin<3, MT=[]; end

% subtract trend if needed
if trend.flag ==1
    G= [ones(size(trend.xy,1), 1), trend.xy(:,1), trend.xy(:,2)];
    beta = G\yVal; 
    dpred = G*beta; 
    yValDetrend= yVal - dpred; 
    trend.Trend_Params=beta;
    trend.Dpred = dpred;
else
    yValDetrend=yVal;
end

% build raw variogram
rawVario = 0.5*L2_distance(yValDetrend,yValDetrend).^2;

% if anisotropy matrix provided, transform to isotropic coordinates
if ~isempty(MT)
   yLocEff=(MT*yLoc')';
else
   yLocEff=yLoc;
end
hEff = L2_distance(yLocEff,yLocEff);

end