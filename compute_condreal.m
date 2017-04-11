function [Dreal] = compute_condreal(xy, XY, d,model, c, lambda, Nreal)
% This function computes conditional realizations of the field v located at
% points xy at locations given in XY, using the data covariance SIG and the
% joint data-estimation covariances sig0, using the lambda values obtained
% from kriging. 

    % parameters
    Nobs = size(xy,1); 
    Nest = size(XY, 1); 
    
    X = [XY; xy]; 
    [COV] = compute_covariance(model, c, X);

    % compute the cholesky factorization of SIG
    [C, p] = chol(COV);
    if p>0, [C,p] = chol(COV+4*sqrt(eps)*eye(size(COV))); end
    if p~=0, error('COV not PD'); end
    
    % Generate conditional realizations
    [Dreal] = deal(zeros(Nest, Nreal)); 
    for k = 1:Nreal
        
        % Generate an unconditional realization
        Vk = C*randn(Nobs+Nest, 1); 
        
        % identify the observation locations
        obs_V = Vk(end-Nobs+1:end,:);
        
        % subtract out the observed from the simulated data
        dV = (d-obs_V);
        
        % calculate the conditional realization
        Dreal(:, k)=Vk(1:Nest)+lambda'*dV;
    end
    
end