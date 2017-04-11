function [C] = compute_covariance(model, c, xy, XY, include_nugget)
% This function returns the covariance matrix given a covariance model with
% parameters c and locations xy. 
% Inputs: 
%   model       - can be Gaussian or exponential
%   c           - model parameters: Generally [sill, range, nugget]
%   xy          - matrix of locations
%   
% Outputs
%   C           - theoretical covariance matrix
%
% Author: Jeremy Maurer, April 17, 2017
% License: MIT

% default parameter values
if nargin < 5, include_nugget=0; end

% compute distance
if numel(xy) == 1
    h = 0; 
elseif nargin < 4
    h = distance_(xy, xy); 
else
    h = distance_(xy, XY); 
end

switch model
    case 'nugget'
        C = zeros(size(h)); 
        C(1) = c; 

    case 'Gaussian'
        % set a lower threshold on covariance and always include nugget
        thresh = (c(1)+c(3))*eps; 
        [~, C] = gaussianVario(c, h, thresh,include_nugget); 
        
    case 'exponential'
        [~, C] = exponentialVario(c, h);

    case 'power'
        [~, C] = linearVario(c, h, include_nugget);

    case 'matern'
        nu = c(3); 
        if include_nugget==1,c = [c(1:2), c(4)]; else, c= c(1:2); end
        [~,C] = maternVario(c, nu, h, include_nugget);
end

end