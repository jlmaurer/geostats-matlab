function [] = plot_variogram(hraw, vraw, N, plotraw)
% This function plots the empirical, binned, and theoretical variograms. 

if nargin < 4
    if nargin < 3
        N = 15; 
    end
    plotraw = 1; 
end

rvar= max(hraw(:))/2;
step = rvar/(N-5); 
xbin1 = [0:step:rvar]; 
xbin2 = linspace(0,ceil(step/3),5); 
xbin = sort([xbin1,xbin2]); 

% compute binned variogram
[hb, vb] = varioexp_(hraw, vraw, xbin);

figure; 
if plotraw==1
    scatter(hraw(:),vraw(:), '.', 'MarkerFaceColor', [.8 .8 .8])
end
hold on
plot(hb,vb, '-*k')
xlabel('Distance')
ylabel('Semi-variance')
title('Raw and Binned variograms')
end

