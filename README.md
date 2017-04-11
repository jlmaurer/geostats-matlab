# geostats-matlab
This repository contains a set of functions and testing scripts for geostatistical analysis, include variogram plotting and 
analysis, kriging, random field simulation, and conditional simulation of data. The main function for analysis is geostatm.m. 
This function accepts a covariance model (can be Gaussian, exponential, spherical, matern, nugget, or power law), as well as 
initial parameter guesses and the data to be interpolated. This package is still under some development, and I have mainly only
tested it on 2-D data sets, with some limited tests in 1D. 

## Dependencies
variogramfit.m and randomfield.m are two optional functions that can be used with geostatm. randomfield was developed by 
Qiqi Wang (qiqi@mit.edu) and Paul G. Constantine (paul.constantine@stanford.edu). It is used to generate random fields with 
a given covariance structure, and can also be used for conditional simulation, although currently I don't use their code to 
do this. Variogramfit is an alternative to lsqcurvefit from Wolfgang Schwanghart which I use by default for estimating the parameters of the theoretical variogram. In my experience, lsqcurvefit typically does better but the difference may be small depending on the situation. Both of these codes can be obtained from the Matlab file Exchange. The only required dependency is 
that testing_geostatm currently calls randomfield. 

## License
All of the codes in this repository except those listed above under dependencies is licensed under the MIT license. It is free to use and distribute. If you use this code for a publication please do include a citation to this repository. Licenses for the dependencies above are included in the repository.

## Example
xy = rand(10, 2); 
data = rand(10,1); 
modelparams.covtype = 'nugget';
modelparams.c0 = [mean(data)];
[X,Y] = meshgrid(linspace(0, 1, 10), linspace(0,1,10)); 
modelparams.XY = [X(:), Y(:)];
modelparams.lb = 0; 
modelparams.ub = 1; 
Nboot = 10; 
Nreal = 10; 
trend = struct('flag', 0);
[Dest, Dsig, Dreal, param, trEst] = geostatm (xy, data, modelparams, trend, Nboot, Nreal)