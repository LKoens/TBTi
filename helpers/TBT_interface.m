function [SO0,SNO] = TBT_interface(epps,rho,rrhop,kappa,r1,r2,r3,t,erho1,erho2,erho3,itau,d, lambda, N, numS, numPhi, tol)
%ACSBT Calculates the results of TBT for a prolate ellipsoid by an
%interface
%NEW APPROACH FOR SYSTEMS REQUIREING HIGH N. AT TESTING STAGE.
%  INPUTS
%  epps     -   slenderness parameter
%  rho(s)   -   cross sectional shape
%  rrhop(s) -   derivative of cross sectional shape times (rho(s)rho'(s))
%  kappa(s) -   curvature of centerline
%  r1(s)    -   centerline in x direction
%  r2(s)    -   centerline in y direction
%  r3(s)    -   centerline in z direction
%  t(s)     -   tangent to centerline
%  erho1(s,phi)-   local radial vector component in x
%  erho1(s,phi)-   local radial vector component in y
%  erho1(s,phi)-   local radial vector component in z
%  d        - distance bellow the interface
%  lambda   - ratio of the viscosities of the two regions
%  itau(s)  -   integrated torsion along the body
%  N        -   Order to expand the solution to
% numS     - number of devisions in s
% numPhi     - number of devisions in Phi
% tol - absolute tolerance for integrals
% 
% OUTPUTS
%  SO0        -   First matrix inversion
%  SNO        -   Subsequent matix inversions
coder.extrinsic('quad2dv')
coder.extrinsic('textprogressbar')

% coefficents for the solutions
lambdaB = (1-lambda)/(1+lambda);
lambda2 = 2 * lambda/(1+lambda);
Amat = [1,0,0;0,1,0;0,0,-1];
Bmat = [lambdaB ,0,0;0,lambdaB,0;0,0,-1];

%surface

S1 = @(s,phi)r1(s)+epps * rho(s) .* erho1(s,phi);
S2 = @(s,phi)r2(s)+epps * rho(s) .* erho2(s,phi);
S3 = @(s,phi)r3(s)+epps * rho(s) .* erho3(s,phi)-d;

ds = 2/(numS); %distance between consecutive s points.
s = linspace(-1+(ds/2),1-(ds/2),numS); %s points to be used for colocation method
dphi = 2 * pi/numPhi;
phi = linspace(-pi+dphi/2,pi-dphi/2,numPhi);
SO0 = zeros(3 * numS * numPhi,3 * numS * numPhi); %solution matrix for next order
SNO = zeros(3 * numS * numPhi,3 * numS * numPhi); %solution matrix for next order

% geometery and velocity
u = zeros(numPhi,numS);
c = zeros(numPhi,numS);
a = zeros(numPhi,numS);
alpha = zeros(numPhi,numS);
tv = zeros(numPhi,numS,3);
for i = 1:numS
    for j = 1:numPhi
        tv(j,i,:) = t(s(i));
        rt = rho(s(i))^2 +sqrt(rho(s(i))^4+4 * rrhop(s(i))^2);
        c(j,i) = sqrt(rt/2);
        u(j,i) =-2 * rrhop(s(i))/rt;
        a(j,i) = 1-epps * rho(s(i)) * kappa(s(i)) * cos(phi(j)-itau(s(i)));
        alpha(j,i) = epps * c(j,i)/a(j,i);
    end
end


%construct zeta_1' and zeta_2' at points s
z1p = Z1p(a,alpha,u);
z2p = Z2p(a,alpha,u);

%construct <1/zeta_1'>  <1/zeta_2'> at points s
z1pib = dphi * sum(1 ./ z1p); 
z2pib = dphi * sum(1 ./ z2p);

% Construct the matrix B.
Bl = zeros(3 * numS);
Bn = zeros(3 * numS);
%textprogressbar('Building B:')
for i = 1 : numS
    Bl(((1:3)+(i-1) * 3),((1:3)+(i-1) * 3)) = t(s(i)).*t(s(i))'/z1pib(i) + (eye(3) - t(s(i)).*t(s(i))')/z2pib(i);
    for j = 1 : numS
        % Perform the integrals of the various contributions of all the
        % kernels. As the domain of integration is shared, and it's 1D, we can
        % compute all of these integrals in one call to integral. Output is
        % returned as a struct with expected field names.
        integrand = @(s2) B_integrand(s2, s(i), r1, r2, r3, rho, epps, d);
        o = unpack_B_integrals(integral(integrand, s(j)-ds/2, s(j)+ds/2, 'ArrayValued', true, 'AbsTol', tol));

        BS = o.Bi * eye(3) +[o.B11 o.B12 o.B13; o.B12 o.B22 o.B23; o.B13 o.B23 o.B33];
        BSm =(o.Biw * eye(3) +[o.B11w o.B12w o.B13w; o.B12w o.B22w o.B23w; o.B13w o.B23w o.B33w]) * Bmat;
        BSDm = -lambda2 * (o.Biwsd * eye(3) -3 * [o.B11wsd o.B12wsd o.B13wsd; o.B12wsd o.B22wsd o.B23wsd; o.B13wsd o.B23wsd o.B33wsd]) * Amat;
        BSSm = -lambda2 * [0,0,o.B13wss;0,0,o.B23wss; -o.B13wss,-o.B23wss,0] * Amat;
        Bn(((1:3)+(i-1) * 3),((1:3)+(j-1) * 3)) = BS + BSm + BSDm + BSSm;
    end
    %textprogressbar(i / numS * 100)
end
%textprogressbar('.')
B = Bl + Bn;


Mi = zeros(3,3,numPhi,numS,numPhi,numS);
for j = 1:numS
    for i = 1:numPhi
        Mi(:,:,i,j,i,j) = (t(s(j)) .* t(s(j))')/z1p(i,j) +((eye(3)-t(s(j)) .* t(s(j))'))/z2p(i,j);
    end
end

MT0 = zeros(3,3,numS,numPhi,numS);
for k = 1:numS
    for l2 = 1:numPhi
        Q1 = reshape(dphi * sum(Mi(:,1,:,:,l2,k),3),3 * numS,1);
        Q2 = reshape(dphi * sum(Mi(:,2,:,:,l2,k),3),3 * numS,1);
        Q3 = reshape(dphi * sum(Mi(:,3,:,:,l2,k),3),3 * numS,1);
        
        A0 = Bl * [Q1, Q2, Q3];

        fb = B \ A0;
        
        res = A0 - Bl * fb;
        MT0(:,1,:,l2,k) = reshape(res(:,1),3,numS);
        MT0(:,2,:,l2,k) = reshape(res(:,2),3,numS);
        MT0(:,3,:,l2,k) = reshape(res(:,3),3,numS);

    end
end

%textprogressbar('Building SO0:')
for k = 1:numS
    for l2 = 1:numPhi
        for j = 1 : numS
            T = t(s(j)) .* t(s(j))';
            for i = 1:numPhi
                MiContrib = (T/z1p(i,j) +((eye(3)-T)/z2p(i,j)));
                SO0((1:3)+3 * (j-1)+3 * numS * (i-1),(1:3)+3 * (k-1)+3 * numS * (l2-1)) = MiContrib*(i==l2 && j==k) - MiContrib * MT0(:,:,j,l2,k);
            end
        end
    end
    %textprogressbar(k / numS * 100)
end
%textprogressbar('.')


%% Next order creation

if N ~= 0
    
    %Create the full matrix for integral equations
    ci  = zeros(numPhi * numS * numPhi * numS,1);
    c11 = zeros(numPhi * numS * numPhi * numS,1);
    c12 = zeros(numPhi * numS * numPhi * numS,1);
    c13 = zeros(numPhi * numS * numPhi * numS,1);
    c22 = zeros(numPhi * numS * numPhi * numS,1);
    c23 = zeros(numPhi * numS * numPhi * numS,1);
    c33 = zeros(numPhi * numS * numPhi * numS,1);
    
    ciw  = zeros(numPhi * numS * numPhi * numS,1);
    c11w = zeros(numPhi * numS * numPhi * numS,1);
    c12w = zeros(numPhi * numS * numPhi * numS,1);
    c13w = zeros(numPhi * numS * numPhi * numS,1);
    c22w = zeros(numPhi * numS * numPhi * numS,1);
    c23w = zeros(numPhi * numS * numPhi * numS,1);
    c33w = zeros(numPhi * numS * numPhi * numS,1);
    
    ciwsd  = zeros(numPhi * numS * numPhi * numS,1);
    c11wsd = zeros(numPhi * numS * numPhi * numS,1);
    c12wsd = zeros(numPhi * numS * numPhi * numS,1);
    c13wsd = zeros(numPhi * numS * numPhi * numS,1);
    c22wsd = zeros(numPhi * numS * numPhi * numS,1);
    c23wsd = zeros(numPhi * numS * numPhi * numS,1);
    c33wsd = zeros(numPhi * numS * numPhi * numS,1);
    
    c13wss = zeros(numPhi * numS * numPhi * numS,1);
    c23wss = zeros(numPhi * numS * numPhi * numS,1);
    
    numEntries = numPhi * numS * numPhi * numS;
    s = linspace(-1+(ds/2),1-(ds/2),numS);
    phi = linspace(-pi+dphi/2,pi-dphi/2,numPhi);
    
    DQ = parallel.pool.DataQueue;
    %textprogressbar('Building SNO:')
    afterEach(DQ, @nUpdateProgressbar);
    p = 1;
    parfor l = 1 : numEntries 
        [k,l2,j,i,sameCell] = getIndex(l,numS,numPhi);
        
        pt = phi(i);
        st = s(j);
        
        a2 = a(i,j);
        t12 = tv(i,j,1);
        t22 = tv(i,j,2);
        t32 = tv(i,j,3);
        r12 = r1(st);
        r22 = r2(st);
        r32 = r3(st);
        c2 = c(i,j);
        u2 = u(i,j); %effective s point of ellipsoid
        r1e = @(s2) a2 .* (s2-u2) .* t12 +r12;
        r2e = @(s2) a2 .* (s2-u2) .* t22 +r22;
        r3e = @(s2) a2 .* (s2-u2) .* t32 +r32;
        rhoe = @(s2) c2 .* sqrt(1-s2.^2);
        S1e = @(s2,phi2) r1e(s2) +epps * rhoe(s2) .* erho1(st,phi2);
        S2e = @(s2,phi2) r2e(s2) +epps * rhoe(s2) .* erho2(st,phi2);
        S3e = @(s2,phi2) r3e(s2) +epps * rhoe(s2) .* erho3(st,phi2);
        s3 = @(s2) u2-st+s2;
        
        % Normal contributions
        C = zeros(7,1);
        if sameCell % Cell with the singularity within it
            
            low = max(s(k)-ds/2,-1+st-u2);
            upp = min(s(k)+ds/2,1+st-u2);
        
            C = C + quad2dv(@(s,phi) c1Integrand(s,phi,st,pt,u2,s3,S1,S2,S3,S1e,S2e,S3e),low,st,phi(l2)-dphi/2,phi(l2)+dphi/2,'AbsTol',tol);
            C = C + quad2dv(@(s,phi) c2Integrand(s,phi,st,pt,u2,s3,S1,S2,S3,S1e,S2e,S3e),st,upp,phi(l2)-dphi/2,phi(l2)+dphi/2,'AbsTol',tol);
            C = C + quad2dv(@(s,phi) c3Integrand(s,phi,st,pt,u2,s3,S1,S2,S3,S1e,S2e,S3e),-1+st-u2,low,-pi,pi,'AbsTol',tol);
            C = C + quad2dv(@(s,phi) c4Integrand(s,phi,st,pt,u2,s3,S1,S2,S3,S1e,S2e,S3e),s(k)-ds/2,low,phi(l2)-dphi/2,phi(l2)+dphi/2,'AbsTol',tol);
            C = C + quad2dv(@(s,phi) c5Integrand(s,phi,st,pt,u2,s3,S1,S2,S3,S1e,S2e,S3e),upp,1+st-u2,-pi,pi,'AbsTol',tol);
            C = C + quad2dv(@(s,phi) c6Integrand(s,phi,st,pt,u2,s3,S1,S2,S3,S1e,S2e,S3e),upp,s(k)+ds/2,phi(l2)-dphi/2,phi(l2)+dphi/2,'AbsTol',tol);
            C = C + quad2dv(@(s,phi) c7Integrand(s,phi,st,pt,u2,s3,S1,S2,S3,S1e,S2e,S3e),low,upp,-pi,phi(l2)-dphi/2,'AbsTol',tol);
            C = C + quad2dv(@(s,phi) c8Integrand(s,phi,st,pt,u2,s3,S1,S2,S3,S1e,S2e,S3e),low,upp,phi(l2)+dphi/2,pi,'AbsTol',tol);

        else % If there is no singularity within the cell we can evaluate the terms directly

            integrand = @(x,y) cStokesletIntegrand(x,y,st,pt,S1,S2,S3,S1e,S2e,S3e);
            C = quad2dv(integrand,s(k)-ds/2,s(k)+ds/2,phi(l2)-dphi/2,phi(l2)+dphi/2,'AbsTol',tol);
            
        end

        ci(l)  = C(1);
        c11(l) = C(2);
        c12(l) = C(3);
        c13(l) = C(4);
        c22(l) = C(5);
        c23(l) = C(6);
        c33(l) = C(7);
        
        % INTERFACE CONTRIBUTIONS
        integrand = @(x,y) CWallIntegrand(x,y,st,pt,S1,S2,S3,S1e,S2e,S3e);
        D = quad2dv(integrand,s(k)-ds/2,s(k)+ds/2,phi(l2)-dphi/2,phi(l2)+dphi/2,'AbsTol',tol); 
        ciw(l)    = D(1);
        c11w(l)   = D(2);
        c12w(l)   = D(3);
        c13w(l)   = D(4);
        c22w(l)   = D(5);
        c23w(l)   = D(6);
        c33w(l)   = D(7);
        ciwsd(l)  = D(8);
        c11wsd(l) = D(9);
        c12wsd(l) = D(10);
        c13wsd(l) = D(11);
        c22wsd(l) = D(12);
        c23wsd(l) = D(13);
        c33wsd(l) = D(14);
        c13wss(l) = D(15);
        c23wss(l) = D(16);
        
        send(DQ, l);
    end
    %textprogressbar('.')

    %% CONSTRUCTING MATRIX SYSTEM.
    C = zeros(3 * numS * numPhi,3 * numS * numPhi); %Construct the matrix C
    
    for l = 1 : numEntries
        [k,l2,j,i,sameCell] = getIndex(l, numS, numPhi);
        
        CS =ci(l) * eye(3) +[c11(l) c12(l) c13(l); c12(l) c22(l) c23(l); c13(l) c23(l) c33(l)];
        CSm =(ciw(l) * eye(3) +[c11w(l) c12w(l) c13w(l); c12w(l) c22w(l) c23w(l); c13w(l) c23w(l) c33w(l)]) * Bmat;
        CSDm = -lambda2 * (ciwsd(l) * eye(3) -3 * [c11wsd(l) c12wsd(l) c13wsd(l); c12wsd(l) c22wsd(l) c23wsd(l); c13wsd(l) c23wsd(l) c33wsd(l)]) * Amat;
        CSSm = -lambda2 * [0,0,c13wss(l);0,0,c23wss(l); -c13wss(l),-c23wss(l),0] * Amat;
        C((1:3)+3 * (j-1)+3 * numS * (i-1),(1:3)+3 * (k-1)+3 * numS * (l2-1)) = CS -Bn(((1:3)+(j-1) * 3),((1:3)+(k-1) * 3)) * dphi+CSm+CSDm+CSSm;
        if sameCell
            C((1:3)+3 * (j-1)+3 * numS * (i-1),(1:3)+3 * (k-1)+3 * numS * (l2-1)) = C((1:3)+3 * (j-1)+3 * numS * (i-1),(1:3)+3 * (k-1)+3 * numS * (l2-1)) -( (t(s(j)) .* t(s(j))') * (z1p(i,j)-Z1(a(i,j),alpha(i,j))) +((eye(3)-t(s(j)) .* t(s(j))')) * (z2p(i,j)-Z2(a(i,j),alpha(i,j))));
        end
    end
    
    
    
    SNO = SO0 * C;
    
end

function nUpdateProgressbar(~)
        %textprogressbar( p / numEntries * 100);
        p = p + 1;
end

end