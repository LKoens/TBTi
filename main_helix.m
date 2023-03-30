addpath(genpath('helpers'))
%% Setup

% Slenderness of the body.
epps = 0.05;

%Helix radius
rh=0.0510937500000000;

% wavenumber
k = 10.5744;

%Helix angle
alpha=sqrt(1-(k*rh)^2);

% Cross-sectional radius; function of arclength.
rho = @(s) (1-s.^(20)).^(0.5);

% rho * d/ds (rho); function of arclength.
rrhop = @(s) -10*s.^(19);

% Centreline curvature; function of arclength.
kappa = @(s)(k^2)*rh;

%theta rotation around the helix axis
theta=pi/2;

% Cartesian components of the centreline; function of arclength.
r1 = @(s) alpha*s;
r2 = @(s)rh*cos(k*s+theta);
r3 = @(s)rh*sin(k*s+theta);

% Tangent vector to the centreline; function of arclength.
t =@(s)[alpha,-k*rh*sin(k*s+theta),k*rh*cos(k*s+theta)] + 0*s;

% Cartesian components of the local radial vector; function of arclength and angle.
erho1 = @(s,phi)k*rh*sin(phi-k*alpha*s);
erho2 = @(s,phi)-cos(k*s+theta).*cos(phi-k*alpha*s)+alpha*sin(k*s+theta).*sin(phi-k*alpha*s);
erho3 = @(s,phi)-sin(k*s+theta).*cos(phi-k*alpha*s)-alpha*cos(k*s+theta).*sin(phi-k*alpha*s);

% Integral of the torsion along the centreline.
itau = @(s)k*alpha*s;

% z offset from the plane boundary.
d = 10000;

% Ratio of viscosity of the two fluid regions.
lambda =0;

% Number of subdivisions of s and Phi.
numS = 30;
numPhi = 30;

% Absolute tolerance for integrals.
tol = 1e-6;

%% Call TBT_interface
tic
[SO0,SNO] = TBT_interface(epps,rho,rrhop,kappa,r1,r2,r3,t,erho1,erho2,erho3,itau,d, lambda, 1, numS, numPhi, tol);
toc
% tic
% [SO0OLD,SNOOLD] = TBT_interfacev3OLD(epps,rho,rrhop,kappa,r1,r2,r3,t,erho1,erho2,erho3,itau,d, lambda, N, numS);
% toc

disp('Iterating')
tic
[R,fs,fsTotal,S,Phi] = Rmat2(SO0,SNO,numS,numPhi,epps,rho,r1,r2,r3,erho1,erho2,erho3);
toc
S = S'; Phi = Phi';

% Plot the final approximation to f in the 6 test cases needed to compute R.
plotF(fsTotal{end},S,Phi)