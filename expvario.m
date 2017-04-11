function [hExp,expVario] = expvario(hEff,rawVario,xBin)
%expvario - calculate the experimental variogram from the raw variogram
%
%SYNTAX: [hExp,expVario] = expvario(hEff,rawVario,xBin)
%
%DESCRIPTION:
%   This function calculates the experimental variogram from the raw variogram
%   using the intervals specified by the vector xBin.  
%
%INPUT PARAMETERS (Required):
%   hEff     - the effective separation distance (incorporating anisotropy) [n^2 x 1]
%   rawVario - raw data semivariances
%   xBin     - vector of the edges of bins.  The values of xBin must be monotonically
%              increasing [a+1 x 1]
%OUTPUT PARAMETERS:
%   hExp     - the mean of the effective distances within the bin [1 x a]
%   expVario - the mean of the raw variograms within the bin  [1 x a]
%EXTERNAL FUNCTIONS USED: distance_
%
%EXAMPLE: 
%   yLoc=20*rand(5,2);
%   yVal=rand(5,1);
%   [hEff,rawVario] = varioraw_(yLoc,yVal)
%   xBin=[0,1,2,5,10];
%   [hExp,expVario] = varioexp_(hEff,rawVario,xBin)
%   plot(hEff,rawVario,'*',hExp,expVario,'k*-');
%
%SEE ALSO: rawvario

%REVISION HISTORY:
%   tae, 3/12/99 created --> original code from the Geostatistics class 
%       taught by Anna Michalak at Stanford University
%   jlm, 4/11/2017 modified, name changed from varioexp_ to expvario

% build the experimental variogram
nBins=length(xBin)-1;
[hExp,expVario] =deal(zeros(1,nBins)); 
for iBin=1:nBins
   iBinContains=find(xBin(iBin)<hEff & hEff<=xBin(iBin+1));
   
   if isempty(iBinContains)
      hExp(iBin)=NaN;
      expVario(iBin)=NaN;
   else
      hExp(iBin) = mean(hEff(iBinContains));
      expVario(iBin) = mean(rawVario(iBinContains));
   end  
end

end