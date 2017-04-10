function [Dest, Dsig, Dreal, param, trEst] = geostatm (xy, data, XY, model,c0,trend, Nest, Nreal)
% This function runs a complete geostatistical analysis workflow for data located at points xy (assuming 
% Cartesian coordinate system) and interpolates onto points XY. The function also can return conditional
% realizations. However, this function won't do all the work for you, because it expects a covariance model,
% initial guesses for the parameters, and whether or not you expect a trend in the data. At this point only 
% linear trends in x and y are included. This function currently only implements ordinary and universal 
% kriging. The function assumes you've already corrected for anisotropy in xy and XY.
%
% Inputs: 
%   xy      - observation locations in 2-D space, assumes a Cartesian coordinate system
%   data    - observation values at each location
%   XY      - estimation points in 2-D space
%   model   - covariance model to use
%   c0      - initial guesses for the parameters to use for the covariance model. 
%   trend   - structure containing the following fields:
%               - flag:     0 or 1, 1 if a trend is to be estimated
%               - xy:       observation locations (only needed if trend.flag ==1). Currenlty only linear trends 
%                           in xy-space are used.
%   Nest    - If provided, Nest is the number of bootstrap estimates of the variogram parameters to perform for 
%               uncertainty analysis.
%   Nreal   - [] (default) if only kriging is desired, otherwise an integer number of conditional realizations 
%               to perform.
%
% Outputs: 
%   Dest    - Kriging estimates at the estimation locations XY
%   Dsig    - Kriging uncertainty at XY
%   Dreal   - Nreal conditional realizations at the estimation locations
%   param   - estimated variogram parameters for the specified covariance model
%   trEst   - estimate parameters for the trend

% specify defaults
if nargin < 8,  Nreal = []; end
if nargin < 7,  Nest = []; end
if nargin < 6,  trend =struct('flag', 0); end

% The function assumes you've already corrected for anisotropy
[h,v, trEst] = rawvario(xy,data(:), [], trend); 
% plot_variogram(h, v)

%% Fitting a variogram
% param = [sill, range, nugget] for Gaussian, exponential
solver = 'lsqcurvefit'; 
% solver = 'variogramfit'; 
if isempty(Nest)
    param = estimate_SVparams(model, h, v, c0);
else
    [~, param] = bootstrap_vario(h, v, model, c0, solver, Nest);
end

% Data Covariance with a specified tolerance
[SIG] = compute_covariance(model, param, xy); 
[sig0] = compute_covariance(model, param, xy, XY); 
[sig2] = compute_covariance(model, param, 0); 

if trend.flag==0
    type = 'ok'; 
elseif trend.flag==1
    type = 'uk'; 
end

%% Kriging
% Currently Ordinary ('ok') and Universal ('uk') kriging 
% are allowed. 

switch type
    case 'uk'
        [Dest, Dsig] = kriging_cpuk(SIG, sig0,data, sig2,xy, XY);
    case 'ok'
        [Dest, Dsig] = kriging_cpok(SIG, sig0,data, sig2);
end

%% Conditional Realizations
if ~isempty(Nreal)
    Dreal = compute_condreal(xy, XY, data, model,param, Nreal); 
else
    Dreal = []; 
end
end
