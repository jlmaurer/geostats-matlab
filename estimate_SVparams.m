function [params,exitflag,output] = estimate_SVparams(model, hraw, vraw, c0, lb, ub)
% This function estimates the best-fitting semivariogram parameters given
% a model and raw variogram. 

% pick the model to use
all_models = {'Gaussian', 'gaussian', 'gauss', 'g', 'exp', 'Exp', 'Exponential', ...
    'exponential', 'e', 'Power', 'power', 'p', 'Spherical', 'spherical',...
    'Nugget', 'nugget', 'n', 'matern'}; 
modvec = strcmp(model, all_models); 
mod = find(modvec==1); 

% specify c0 if none supplied
nparams = [3*ones(1,length(all_models)-1), 1, 1, 1];
if nargin < 4, c0 = zeros(1,nparams(mod)); end

if mod == 1 || mod == 2 || mod == 3 || mod == 4
    f = @gaussianVario; 

elseif mod == 5 || mod ==6 || mod ==7 || mod == 8 || mod ==9
    % Exponential
    if numel(c0)==3
        f = @(param,h) exponentialVario(param,h, 1);    
    else
        f = @exponentialVario;
    end
    
elseif mod == 10 || mod ==11 || mod ==12 
    % Power Law
    if numel(c0)==3
        f = @(c, h) c(1)*(h.^c(2)) + c(3)*(h~=0);
    else
        f = @(c, h) c(1)*(h.^c(2)); 
    end
    
elseif mod == 13 || mod ==14
    % Spherical
    if numel(c0)==3
        f = @(c, h) c(1)*(1.5*h./c(2) - 0.5*(h.^3./c(2).^3)) + c(3)*(h~=0);
    else
        f = @(c, h) c(1)*(1.5*h./c(2) - 0.5*(h.^3./c(2).^3));
    end
    
elseif mod == 15 || mod ==16 || mod==17
    % Nugget
    params = mean(vraw); 
    exitflag = []; 
    output = []; 
    return
    
elseif mod==18
    nu = c0(3); 
    if numel(c0)==4
        c0 = [c0(1:2), c0(4)]; 
        f = @(c, h) maternVario(c,nu, h, 1); 
    else
        c0 = [c0(1:2)]; 
        f = @(c, h) maternVario(c, nu, h);
    end

else
    error('Please specify a model: Gaussian, Exponential, Power, Nugget, or Spherical')
end

[params,~,~,exitflag,output] = lsqcurvefit(f,c0,hraw,vraw,lb,ub);

if mod==18, params = [params(1:2), nu, params(3)]; end
end
