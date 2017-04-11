function [Dest, Dsig, Dreal, param, trEst] = geostatm (xy, data, modelparams, trend, Nboot, Nreal)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%   model   - covariance model to use. 'Gaussian,' 'Exponential,' 'nugget,'
%               'linear,' 'matern.'
%   c0      - initial guesses for the parameters to use for the covariance model. 
%   trend   - structure containing the following fields:
%               - flag:     0 or 1, 1 if a trend is to be estimated
%               - xy:       observation locations (only needed if trend.flag ==1). Currenlty only linear trends 
%                           in xy-space are used.
%   Nboot   - If provided, Nest is the number of bootstrap estimates of the variogram parameters to perform for 
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
%
% Author: Jeremy  Maurer, copyright 2017. MIT License. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% specify defaults
if nargin < 6,  Nreal = []; end
if nargin < 5,  Nboot = []; end
if nargin < 4,  trend =struct('flag', 0); end

XY = modelparams.XY; 
model = modelparams.covtype; 
c0 = modelparams.c0; 
lb = modelparams.lb; 
ub = modelparams.ub; 

% The function assumes you've already corrected for anisotropy
[h,v, trEst] = rawvario(xy,data(:), [], trend); 

% uncomment below to plot (empirical) variogram. Helpful for
% troubleshooting.
% plot_variogram(h, v)

%% Fitting a variogram
% param = [sill, range, nugget] for Gaussian, exponential
solver = 'lsqcurvefit'; 
% solver = 'variogramfit'; 
if isempty(Nboot)
    if strcmp(solver, 'lsqcurvefit')==1
            param = estimate_SVparams(model, h, v, c0, lb, ub);
    elseif strcmp(solver, 'variogramfit')==1
        nlag = 20; 
        [param(2),param(1),param(3)] = variogramfit(h(:),v(:),c0(2),c0(1),nlag,...
            'model', 'gaussian', 'nugget', c0(3));
    else
        error('Please use a valid solver');
    end

else
    [~, param] = bootstrap_vario(h, v, model, c0, solver, lb, ub, Nboot);
end

% Data Covariance with a specified tolerance
[SIG] = compute_covariance(model, param, xy, xy, 1); 
[sig0] = compute_covariance(model, param, xy, XY,0); 
[sig2] = compute_covariance(model, param, 0, 0, 1); 

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
        [Dest, Dsig, lambda] = kriging_cpuk(SIG, sig0,data, sig2,xy, XY);
    case 'ok'
        [Dest, Dsig, lambda] = kriging_cpok(SIG, sig0,data, sig2);
end

%% Conditional Realizations
if ~isempty(Nreal)
    Dreal = compute_condreal(xy, XY, data, model,param, lambda, Nreal); 
else
    Dreal = []; 
end
end
