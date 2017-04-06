function [Dest, Dsig, Dreal, param, trEst] = geostatm (xy, data, XY, model,c0,trend, type, Nreal)
% This function runs a complete geostatistical analysis workflow for data located at points xy (assuming 
% Cartesian coordinate system) and interpolates onto points XY. The function also can return conditional
% realizations. However, this function won't do all the work for you, because it expects a covariance model,
% initial guesses for the parameters, and whether or not you expect a trend in the data. At this point only 
% linear trends in x and y are included. This function currently only implements ordinary and universal 
% kriging. 
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
%   type    - type of kriging, can be 'ok' for ordinary kriging or 'uk' for universal kriging (kriging with a trend)
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
if nargin < 7, type = 'ok'; end
if nargin < 6, trend =struct('flag', 0); end

[h,v, trEst] = varioraw_(xy,data(:), [], trend); 
% plot_variogram(h, v)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Two choices for fitting a variogram
% params: sill, range, nugget
nboot = 0; 
lbub = [0, max(v)/2; 0, max(h)/2; 0, max(max(v)/4, 2)];
[~, best_p1] = bootstrap_vario(h, v, model, c0, 1, lbub, nboot);
% [~, best_p2] = bootstrap_vario(h, v, model, c0, 2, lbub, nboot);
param = best_p1;
% graph_correlations([p1.s, p1.n, p1.r], 2, {'Sill', 'Nugget', 'Range'}, 0, 0); 
% graph_correlations([p2.s, p2.n, p2.r], 2, {'Sill', 'Nugget', 'Range'}, 0, 0); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hl = linspace(0, 100); 
[V, C] = Gaussian_vario(param, hl);

% Data Covariance with a specified tolerance
tol = 1e-5*max(C); 
[~,SIG] = Gaussian_vario(param, h, tol); 

h0 = distance_(xy, XY); 
[~,sig0] = Gaussian_vario(param, h0, tol, 0); 
[~,sig2] = Gaussian_vario(param, 0, tol, 1); 

switch type
    case 'uk'
        [Dest, Dsig, lambda] = kriging_cpuk(SIG, sig0,xy, XY, data, sig2);
    case 'ok'
        [Dest, Dsig, lambda] = kriging_cpok(SIG, sig0,data, sig2);
    case 'sk'
        %TODO: need to implement simple kriging
        [Dest, Dsig, lambda] = kriging_cpsk(SIG, sig0,data, sig2);
end

if ~isempty(Nreal)
    Dreal = compute_condreal(SIG, sig0, data, lambda, Nreal); 
else
    Dreal = []; 
end
end

function [hEff,rawVario, trend] = varioraw_(yLoc,yVal, MT,trend)

%REVISION HISTORY:
%   tae, 3/10/99 created. See stand-alone function

% assign empty matrices to unspecified input parameters
if nargin<4,trend = struct('flag',0); end
if nargin<3, MT=[]; end
% if nargin<3, yErr=[]; end

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
rawVario = 0.5*distance_(yValDetrend,yValDetrend).^2;
if ~isempty(MT)
   yLocEff=(MT*yLoc')';              % tranform to generalized distance
else
   yLocEff=yLoc;
end
hEff = distance_(yLocEff,yLocEff);

end

function retVal=VerifyParams_(yLoc,yVal)
% (description of error checking)

[n1,d2]=size(yLoc);
[n2,a1]=size(yVal);
   
% check that the number of measurements are consistant
if n1~=n2
   warning('inconsistant number of measurements: check yLoc and yVal');   
   retVal=0; return;
end;
   
retVal=1; % no errors found
return;
end

function [Dest, Dsig, lambda] = kriging_cpok(SIG, sig0, v, sig2)
% This function returns weights to use for continous-part ordinary Kriging, 
% using the covariance matrices SIG (covariance of the
% observation locations) and sig0 (covariance of the observation locations
% with the measurement locations).All of the estimation locations are
% computed simultaneously. 

    Nobs = size(SIG,2); 
    Nest = size(sig0, 2); 

    A= [SIG, ones(Nobs,1); ones(1,Nobs), 0]; 
    B = [sig0; ones(1,Nest)]; 

    lambdanu = A\B; 
    lambda = lambdanu(1:end-1,:); 
    nu = -lambdanu(end,:); 
    
    % compute predicted velocity and uncertainty
    Dest = lambda'*v; 
    Dsig = sqrt(sig2 - diag(lambda'*sig0) + nu(:)); 
end

function [Dest, Dsig, lambda] = kriging_cpuk(SIG, sig0,xy, XY, v, sig2, dim)
% This function returns weights to use for continous-part Kriging
% w/variable mean, using the covariance matrices SIG (covariance of the
% observation locations) and sig0 (covariance of the observation locations
% with the measurement locations). All of the estimation locations XY (a 
% Nestx2 matrix) are computed simultaneously. xy is the observed locations
% (should be Nobs x 2). 

    Nobs = size(SIG,2); 
    Nest = size(sig0, 2); 
    
    if nargin < 7
        dim=2;
    end
    
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
    
    A= [SIG, F; F', zeros(size(F,2))]; 
    B = [sig0; f0]; 

    lambdanu = A\B; 
    lambda = lambdanu(1:size(sig0,1),:); 
    nu = -lambdanu(size(sig0,1)+1:end,:); 
    fnu = sum(f0.*nu, 1); 
    
    Dest = lambda'*v; 
    Dsig = sqrt(sig2 - diag(lambda'*sig0) + fnu(:)); 
end

function [Dreal] = compute_condreal(SIG, sig0, v, lambda, Nreal)
% This function computes conditional realizations of the field v located at
% points xy at locations given in XY, using the data covariance SIG and the
% joint data-estimation covariances sig0, using the lambda values obtained
% from kriging. 
    Nobs = size(SIG,2); 
    Nest = size(sig0, 2); 

    C = chol(SIG); 
    [Dreal] = deal(zeros(Nest, Nreal)); 
    for k = 1:Nreal
        Vk = C*randn(Nobs+Nest, 1); 
        obs_V = Vk(end-Nobs+1:end,:);
        dV = (v-obs_V);
        for loop=1:Nest
            Dreal(loop, k)=Vk(loop)+lambda(:,loop)'*dV;
        end
    end
    
end

function [params,exitflag,output] = estimate_SVparams(model, hraw, vraw, hcutoff, c0)
% This function estimates the best-fitting semivariogram parameters given
% a model and raw variogram. 

%%%%%%%%%%%%%%
%%TODO: make the nugget optional addition to other models
%%%%%%%%%%%%%%

% pick the model to use
all_models = {'Gaussian', 'gaussian', 'gauss', 'g', 'exp', 'Exp', 'Exponential', ...
    'exponential', 'e', 'Power', 'power', 'p', 'Spherical', 'spherical',...
    'Nugget', 'nugget', 'n'}; 
modvec = strcmp(model, all_models); 
mod = find(modvec==1); 

% specify c0 if none supplied
nparams = [3*ones(1,length(all_models)-1), 1, 1, 1];
if nargin < 5
    c0 = zeros(1,nparams(mod));
end

% specify search bounds
lb = [0, 0, 0]; 
ub = [max(vraw), max(hraw)/2, max(vraw)]; 
ind = hraw < hcutoff; 

if mod == 1 || mod == 2 || mod == 3 || mod == 4
    f = @Gaussian_vario; 
    [params, ~, ~, exitflag, output] = lsqcurvefit(f, c0, hraw(ind), vraw(ind), lb, ub); 
    
elseif mod == 5 || mod ==6 || mod ==7 || mod == 8 || mod ==9
    % Exponential
    f = @Exponential_vario;
    [params, ~, ~, exitflag, output] = lsqcurvefit(f, c0, hraw(ind), vraw(ind), lb, ub); 
    
elseif mod == 10 || mod ==11 || mod ==12 
    % Power Law
    f = @(c, h) c(1)*(h.^c(2)) + c(3)*(h~=0);
    [params, ~, ~, exitflag, output] = lsqcurvefit(f, c0, hraw(ind), vraw(ind), lb, ub); 
    
elseif mod == 13 || mod ==14
    % Spherical
    f = @(c, h) c(1)*(1.5*h./c(2) - 0.5*(h.^3./c(2).^3)) + c(3)*(h~=0);
    [params, ~, ~, exitflag, output] = lsqcurvefit(f, c0, hraw(ind), vraw(ind), lb, ub); 
    
elseif mod == 15 || mod ==16
    % Nugget
    params = mean(vraw); 
    exitflag = []; 
    output = []; 

else
    error('Please specify a model: Gaussian, Exponential, Power, Nugget, or Spherical')
end
    
end

function [var, cov] = Gaussian_vario(param, h, thresh,include_nugget)
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

function h = distance_(x,y)

%REVISION HISTORY:
%   pkk, 5/13/97
%   tae, 6/26/98 (minor changes)
%   jlm, 3/16/2017

[n1,D] = size(x);
[n2,D2] = size(y);

if D~=D2
   error('ERROR in DISTANCE_: locations must have same number of dimensions (columns)')
end
h = zeros(n1,n2);
for id = 1:D
   h = h + (x(:,id)*ones(1, n2)-ones(n1,1)*y(:,id)').^2;
end
h = sqrt(h);
end