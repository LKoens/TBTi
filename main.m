addpath(genpath('helpers'))
%% Setup

% Slenderness of the body.
epps = 1/10;

% Cross-sectional radius; function of arclength.
rho = @(s) (1-s.^2).^(0.5);

% rho * d/ds (rho); function of arclength.
rrhop = @(s) -s;

% Centreline curvature; function of arclength.
kappa = @(s) 0*s;

% Cartesian components of the centreline; function of arclength.
r1 = @(s) s;
r2 = @(s) 0*s;
r3 = @(s) 0*s;

% Tangent vector to the centreline; function of arclength.
t = @(s) [1;0;0] + 0*s;

% Cartesian components of the local radial vector; function of arclength and angle.
erho1 = @(s,phi) 0*s;
erho2 = @(s,phi) cos(phi);
erho3 = @(s,phi) sin(phi);

% Integral of the torsion along the centreline.
itau = @(s) 0*s;

% z offset from the plane boundary.
d = (1+1)*epps;

% Ratio of viscosity of the two fluid regions.
lambda = 10000;

% Number of subdivisions of s and Phi.
numS = 10;
numPhi = 200;

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