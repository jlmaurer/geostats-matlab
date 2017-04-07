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