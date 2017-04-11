function [p,best_param] = bootstrap_vario(h, v, model, c0, solver, lb, ub, Nboot)
% This function bootstraps a variogram solver to get an idea of the
% uncertainty in the estimate. 

nlag = 20; 

h = h(:); 
v = v(:); 
nparam = length(c0); 

param = zeros(Nboot, nparam); 

for loop = 1:Nboot
    [hboot,index] = datasample(h, length(h)); 
    vboot = v(index);
    
    % run the solver
    if strcmp(solver, 'lsqcurvefit')==1
        param(loop,:) = estimate_SVparams(model, hboot, vboot, c0, lb, ub);
    elseif strcmp(solver, 'variogramfit')==1
        if nparam == 3
            [param(loop,2),param(loop, 1),param(loop,3)] = variogramfit(hboot,...
                vboot,c0(2),c0(1),nlag,'model', model, 'nugget', c0(3));
        elseif nparam==2
            [param(loop,2),param(loop, 1),param(loop,3)] = variogramfit(hboot,...
                vboot,c0(2),c0(1),nlag,'model', model);
        elseif nparam==1
            param(loop) = mean(vboot); 
        else
            error('This condition is not yet handled using variogramfit')
        end
    else
        error('Invalid solver')
    end
end

if nparam==3
    p.s = param(:,1); 
    p.r = param(:,2); 
    p.n = param(:,3); 
    names = {'Sill', 'Range', 'Nugget'}; 
elseif nparam==2
    p.s = param(:,1); 
    p.r = param(:,2); 
    names = {'Sill', 'Range'}; 
elseif nparam==1
    p.n = param(:); 
    names = {'Nugget'}; 
else
    error('This condition is not yet handled')
end

if strcmp(solver, 'lsqcurvefit')==1
    best_param = estimate_SVparams(model, h, v, c0, lb, ub);
elseif strcmp(solver, 'variogramfit')==1
    [best_param(2),best_param(1),best_param(3)] = variogramfit(h,v,c0(2),c0(1),nlag,'model', 'gaussian', 'nugget', c0(3));
else
    error('Invalid solver id')
end

graph_correlations(param, 2, names, 0, 0); 

end