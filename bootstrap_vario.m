function [p,best_param] = bootstrap_vario(h, v, model, c0, solver, Nboot)
% This function bootstraps a variogram solver to get an idea of the
% uncertainty in the estimate. 
if nargin < 7, Nboot = 1000;  end
nlag = 20; 

h = h(:); 
v = v(:); 
[all_sills, all_ranges, all_nuggets] = deal(zeros(Nboot,1)); 
for loop = 1:Nboot
    [hboot,index] = datasample(h, length(h)); 
    vboot = v(index);
    
    % run the solver
    param= zeros(1,3);
    if strcmp(solver, 'lsqcurvefit')==1
        param = estimate_SVparams(model, hboot, vboot, c0, max(h(:))/4);
    elseif strcmp(solver, 'variogramfit')==1
        [param(2),param(1),param(3)] = variogramfit(hboot,vboot,c0(2),c0(1),nlag,'model', 'gaussian', 'nugget', c0(3));
    else
        error('Invalid solver')
    end
    all_sills(loop) = param(1); 
    all_ranges(loop) = param(2); 
    all_nuggets(loop) = param(3); 
end
p.s = all_sills; 
p.r = all_ranges; 
p.n = all_nuggets; 

if strcmp(solver, 'lsqcurvefit')==1
    best_param = estimate_SVparams(model, h, v, c0, max(h(:))/4);
elseif strcmp(solver, 'variogramfit')==1
    [best_param(2),best_param(1),best_param(3)] = variogramfit(h,v,c0(2),c0(1),nlag,'model', 'gaussian', 'nugget', c0(3));
else
    error('Invalid solver id')
end

graph_correlations([p.s, p.n, p.r], 2, {'Sill', 'Nugget', 'Range'}, 0, 0); 

end