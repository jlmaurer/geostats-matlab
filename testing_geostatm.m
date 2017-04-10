%% test geostats function on simple field
clear

%% Create synthetic data and plot

% parameters
sd = 1e-1; 
ptru = [3, 10, sd];     % sill, range, nugget
s1 = 0; 
s2 = 0; 

maxX = 50; 
maxY = 50; 

% observation locations
nxy = 100; 
xy = [maxX*rand(nxy,1), maxY*rand(nxy, 1)]; 
% h = distance_(xy, xy); 

[X,Y] = meshgrid(linspace(0, maxX, maxX+1), linspace(0, maxY, maxY+1));
XY = [X(:), Y(:)]; 
allx = [XY; xy]; 
H1 = distance_(allx, allx); 
% H0 = distance_(XY,XY); 
% nXY = (maxX+1)*(maxY+1); 

% compute true data on grid
[~,c] = gaussianVario(ptru, H1); 
corr = struct('C', c, 'name', 'gauss', 'c0', ptru(2), 'sigma', ptru(1));
[F,KL] = randomfield(corr,allx);
dtru = F + s1*allx(:,1) + s2*allx(:,2); 
dtest = dtru(1:length(XY)); 
dobs = dtru(length(XY)+1:end); 
dnoisy = dobs + ptru(3)*randn(size(dobs));
sigscale = (dnoisy-dobs)./ptru(3);
% 
% figure; 
% scatter(XY(:,1), XY(:,2), [], dtest, 'filled')
% hold on
% scatter(xy(:,1), xy(:,2), [], dnoisy, 'filled')

%% Geostatistical Inference
c0 = [2, 15, 0.01]; 

% now call geostatm
trend = struct('flag', 0); 
[~,ind] = sort(xy(:,1));
xy = xy(ind,:); 
dobs = dobs(ind); 
[Vest, Vsig, ~, param, trEst] = geostatm(xy, dobs, XY, 'Gaussian', c0, trend, 2);

Ve = reshape(Vest, maxX+1, maxY+1);
Vs = reshape(Vsig, maxX+1, maxY+1);

Residuals = reshape(Vest - dtest, maxX+1, maxY+1); 
Rn=abs(.5*Residuals./Vs);
[~,sig2] = gaussianVario(param, 0, [], 1); 
totalMisfit = (1/sig2)*Residuals(:)'*eye(length(XY))*Residuals(:);  

figure;
imagesc([0, maxX],[0, maxY],Rn)
caxis([0,1])
colorbar
hold on
scatter(xy(:,1),xy(:,2),[],abs(sigscale),'filled')
set(gca,'ydir','normal')

disp('True parameter values are:')
ptru
disp('Estimated are:')
param

%% Now test other methods

% nearest-neighber default
Vqn = griddata(xy(:,1), xy(:,2), dobs, XY(:,1), XY(:,2), 'nearest'); 
tMn = (1/sig2)*(Vqn-dtest)'*eye(length(XY))*(Vqn-dtest);  

% linear
Vql = griddata(xy(:,1), xy(:,2), dobs, XY(:,1), XY(:,2), 'linear'); 
ind = isnan(Vql); 
Vql(ind) = Vqn(ind); 
tMl = (1/sig2)*(Vql-dtest)'*eye(length(XY))*(Vql-dtest);  

% cubic spline
Vqs = griddata(xy(:,1), xy(:,2), dobs, XY(:,1), XY(:,2), 'cubic');
ind = isnan(Vqs); 
Vqs(ind) = Vqn(ind); 
tMs = (1/sig2)*(Vqs-dtest)'*eye(length(XY))*(Vqs-dtest);  

% natural
Vqnat = griddata(xy(:,1), xy(:,2), dobs, XY(:,1), XY(:,2), 'natural');
ind = isnan(Vqnat); 
Vqnat(ind) = Vqn(ind); 
tMnat = (1/sig2)*(Vqnat-dtest)'*eye(length(XY))*(Vqnat-dtest);  

% v4
Vq4 = griddata(xy(:,1), xy(:,2), dobs, XY(:,1), XY(:,2), 'v4');
ind = isnan(Vq4); 
Vq4(ind) = Vqn(ind); 
tM4 = (1/sig2)*(Vq4-dtest)'*eye(length(XY))*(Vq4-dtest);  

% thin-plate splines
[st,p] = tpaps(xy', dobs');
avals = fnval(st,XY');
tMtp = (1/sig2)*(avals'-dtest)'*eye(length(XY))*(avals'-dtest);  

%% Show results
disp(['Total misfit for Kriging is ', num2str(totalMisfit)])
disp(['Total misfit for nearest-neighber is ',num2str(tMn)])
disp(['Total misfit for linear interp is ',num2str(tMl)])
disp(['Total misfit for cubic interp is ',num2str(tMs)])
disp(['Total misfit for natural interp is ',num2str(tMnat)])
disp(['Total misfit for v4 method is ',num2str(tM4)])
disp(['Total misfit for thin-plate splines is ',num2str(tMtp)])
