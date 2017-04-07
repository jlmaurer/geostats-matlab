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
if nargin < 5, include_nugget=1; end

% compute distance
if numel(xy) == 1
    h = 0; 
elseif nargin < 4
    h = L2_distance(xy, xy); 
else
    h = L2_distance(xy, XY); 
end

switch model
    case 'nugget'
        C = zeros(size(h)); 
        C(1) = c; 

    case 'Gaussian'
        
        % set a lower threshold on covariance and always include nugget
        thresh = (c(1)+c(3))*eps; 
        [~, C] = Gaussian_vario(c, h, thresh,1); 
        
    case 'exponential'
        [~, C] = Exponential_vario(c, h);

    case 'linear'
        [~, C] = linearVario(c, h, include_nugget);

    case 'matern'
        [~,C] = maternVario(c, h, include_nugget);
end

end