function [Mt] = compute_anisotropy_matrix (A, angleinrad)
% This function computes a correction matrix for anisotropy to be used with
% the varioraw_ function. Angle is specified in radians and stretch factor
% A is unitless.

Mt = [1, 0; 0, A]*[cos(angleinrad), sin(angleinrad); -sin(angleinrad), cos(angleinrad)]; 
end
