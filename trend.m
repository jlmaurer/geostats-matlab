function varargout = trend (varargin)
%% Test FMB, PMB, and regular MCMC for simple 1-D,
% multiple station 1-D.

[varargout{1:nargout}]=feval(varargin{:}); 
end

function [vnew, params] = remove_trend (xy, v, flag)
% This function handles the non-linear arctangent trend 
% that is present in cross-fault velocities. 
%
% Input: 
%   xy          - observation locations
%   v           - observation values (Mx2 matrix of velocities)

    vd = var(v); 
    if vd(1) > vd(2)
        component=1;
        other = 2;
    else
        component = 2; 
        other = 1; 
    end

    vtrend = v(:,component);
    trend_dir = 2; 
    xytrend = xy(:,trend_dir); 

    switch flag
        case 'subtract'
            [v_detrend, params] = trend('subtract',xytrend, vtrend);

        case 'transform'
            [v_detrend, params] = trend('transform', xytrend, vtrend);

        case 'nothing'
            v_detrend = vtrend; 
            params = []; 

        case 'linear'
            [vfit, params]= trend('fitpwlin', xytrend, vtrend);  
            v_detrend = vtrend - vfit; 

    end

    vnew = [v_detrend,v(:,other)];
end

function [v_out, params] = subtract (xytrend, vtrend, params)
% This function fits an arctangent to the data and subtracts it out if
% params is not specified. If params is specified, the function adds an
% arctangent trend with the specified parameters back to the data and
% returns the result. 

    fyatan = @(c, x) c(1) - c(2)*atan((x - c(3))./c(4));
    if nargin < 3
        c0 = [0, range(vtrend)/2, 0, range(xytrend)/2]; 
        w = abs(xytrend-mean(xytrend)); w = w./sum(w); 

        % use statistics toolbox
        nlm = fitnlm(xytrend, vtrend,fyatan,c0, 'Weight', w);
        params = nlm.Coefficients.Estimate;
        v_out = vtrend - predict(nlm, xytrend); 
        
    else
        % add trend back to data
        v_out = vtrend + fyatan(params, xytrend); 
    end
end

function [v_out, params] = transform (xytrend, vtrend, params)
% This function returns an tanget transformation of the data and if
% params is not specified. If params is specified, the function takes an
% arctangent transformation with the specified parameters and returns the 
% result. 
    fyatan = @(c, x) c(1) - c(2)*atan((x - c(3))./c(4));
    if nargin < 3
        c0 = [0, range(vtrend)/2, 0, range(xytrend)/2]; 
        w = abs(xytrend-mean(xytrend)); w = w./sum(w); 
%         w = ones(size(xytrend)); w = w./sum(w); 

        % use statistics toolbox
        nlm = fitnlm(xytrend, vtrend,fyatan,c0, 'Weight', w);
        params = nlm.Coefficients.Estimate;
        v_out = tan((vtrend - params(1))./params(2)); 
        
    else
        % add trend back to data
        v_out = atan(vtrend)*params(2) + params(1);  
    end
end
