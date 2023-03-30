function [Q,ERRBND] = quad2dv(FUN,A,B,c,d,varargin)
%QUAD2D    Numerically evaluate double integral over a planar region. Q =
%   QUAD2DV(FUN,A,B,C,D) approximates the integral of FUN(X,Y) over the planar
%   region A <= X <= B and C(X) <= Y <= D(X). FUN is a function handle, C and
%   D may each be a scalar or a function handle. FUN returns pages of output.

%   All input functions must be vectorized:
%   The function Z=FUN(X,Y) must accept 2D matrices X and Y of the same
%   size and return a matrix Z of corresponding values. The functions
%   YMIN=C(X) and YMAX=D(X) must accept matrices and return matrices of the
%   same size with corresponding values.
%
%   [Q,ERRBND] = QUAD2D(...). ERRBND is an approximate upper bound on the
%   absolute error, |Q - I|, where I denotes the exact value of the
%   integral.
%
%   [Q,ERRBND] = QUAD2D(FUN,A,B,C,D,PARAM1,VAL1,PARAM2,VAL2,...) performs
%   the integration as above with specified values of optional parameters:
%
%   'AbsTol', absolute error tolerance
%   'RelTol', relative error tolerance
%
%       QUAD2D attempts to satisfy ERRBND <= max(AbsTol,RelTol*|Q|). This
%       is absolute error control when |Q| is sufficiently small and
%       relative error control when |Q| is larger. A default tolerance
%       value is used when a tolerance is not specified. The default value
%       of 'AbsTol' is 1e-5. The default value of 'RelTol' is
%       100*eps(class(Q)). This is also the minimum value of 'RelTol'.
%       Smaller 'RelTol' values are automatically increased to the default
%       value.
%
%   'MaxFunEvals', maximum number of evaluations of FUN allowed
%
%       The 'MaxFunEvals' parameter limits the number of vectorized calls
%       to FUN. The default is 2000.
%
%   'FailurePlot', generate a plot if MaxFunEvals is reached
%
%       Setting 'FailurePlot' to TRUE generates a graphical representation
%       of the regions needing further refinement when MaxFunEvals is
%       reached. No plot is generated if the integration succeeds before
%       reaching MaxFunEvals. The default is FALSE.
%
%   'Singular', problem may have boundary singularities
%
%       With 'Singular' set to TRUE, QUAD2D will employ transformations to
%       weaken boundary singularities for better performance. The default
%       is TRUE.
%
%   Consider using INTEGRAL2 instead of QUAD2D. INTEGRAL2 is similar to
%   QUAD2D but also supports improper integrals.
%
%   Example:
%   Integrate y*sin(x)+x*cos(y) over pi <= x <= 2*pi, 0 <= y <= pi.
%   The true value of the integral is -pi^2.
%
%      Q = quad2d(@(x,y) y.*sin(x)+x.*cos(y),pi,2*pi,0,pi)
%
%   Example:
%   Integrate 1./(sqrt(x+y).*(1+x+y).^2 over the triangle 0 <= x <= 1,
%   0 <= y <= 1-x. The integrand is infinite at (0,0). The true value of
%   the integral is pi/4 - 1/2.
%
%       fun = @(x,y) 1./( sqrt(x + y) .* (1 + x + y).^2 )
%
%   % In Cartesian coordinates:
%
%       ymax = @(x) 1 - x
%       Q = quad2d(fun,0,1,0,ymax)
%
%   % In polar coordinates:
%
%       polarfun = @(theta,r) fun(r.*cos(theta),r.*sin(theta)).*r
%       rmax = @(theta) 1./(sin(theta) + cos(theta))
%       Q = quad2d(polarfun,0,pi/2,0,rmax)
%
%   Class support for inputs A, B, C, D, and the output of FUN:
%      float: double, single
%
%   See also INTEGRAL2, INTEGRAL, INTEGRAL3, DBLQUAD, TRIPLEQUAD, QUADGK,
%   TRAPZ, FUNCTION_HANDLE, ARRAYFUN.

%   Based on "TwoD" by Lawrence F. Shampine.
%   Ref: L.F. Shampine, "Matlab Program for Quadrature in 2D",
%   Appl. Math. Comp., 202 (2008) 266-274.

%   Copyright 2008-2021 The MathWorks, Inc.

%   Variables accessed in nested functions are written in all caps.
narginchk(5,inf);
if isa(c,'function_handle')
    phiBvar = c;
elseif isfloat(c) && isscalar(c) && isfinite(c)
    phiBvar = @(x) c*ones(size(x));
    end
if isa(d,'function_handle')
    phiTvar = d;
elseif isfloat(d) && isscalar(d) && isfinite(d)
    phiTvar = @(x) d*ones(size(x));
end
USER_SUPPLIED_RELTOL = false;
p = inputParser;
p.addParameter('AbsTol',1e-5,@validateAbsTol);
p.addParameter('RelTol',0,@validateRelTol);
p.addParameter('Singular',true,@validateSingular);
p.addParameter('MaxFunEvals',2000,@validateMaxFunEvals);
p.addParameter('FailurePlot',false,@validateFailurePlot);
p.parse(varargin{:});
optionStruct = p.Results;
ATOL = optionStruct.AbsTol;
RTOL = optionStruct.RelTol;
SINGULAR = logical(optionStruct.Singular);
MaxFunEvals = optionStruct.MaxFunEvals;
if SINGULAR
    thetaL = 0;
    thetaR = pi;
    phiB = 0;
    phiT = pi;
else
    thetaL = A;
    thetaR = B;
    phiB = 0;
    phiT = 1;
end
AREA = (thetaR - thetaL)*(phiT - phiB);
% Gauss-Kronrod (3,7) pair with degrees of precision 5 and 11.
NODES = [ -0.9604912687080202, -0.7745966692414834, -0.4342437493468026, ...
    0, 0.4342437493468026, 0.7745966692414834, 0.9604912687080202 ];
NNODES = length(NODES);
ONEVEC = ones(2*NNODES,1);
NARRAY = [NODES+1,NODES+3]/4;
WT3 = [ 0, 5/9, 0, 8/9, 0, 5/9, 0];
WT7 = [ 0.1046562260264672, 0.2684880898683334, 0.4013974147759622, ...
    0.4509165386584744, 0.4013974147759622, 0.2684880898683334, ...
    0.1046562260264672 ];
VTSTIDX = [16,74,132;27,81,124]; % Some indices between 1 and 4*NNODES^2.
FIRSTFUNEVAL = true;
NFE = 0;
% Compute initial approximations on four subrectangles. Initialize RECTLIST
% of information about subrectangles for which the approximations are not
% sufficiently accurate. NRECTS is the number of subrectangles that remain
% to be processed. ERRBND is a bound on the error.
[Qsub,esub] = tensor(thetaL,thetaR,phiB,phiT);
Q = sum(Qsub,2);
NUMOUTPUTS = size(Qsub,1);
% Now that we know the output class, check RTOL.
outcls = class(Q);
EPS100 = 100*eps(outcls);
if RTOL < EPS100
    RTOL = EPS100;
end
if isa(Q,'double')
    % Single RTOL or ATOL should not force any single precision
    % computations.
    RTOL = double(RTOL);
    ATOL = double(ATOL);
end
rtold8 = max(RTOL/8,EPS100);
atold8 = ATOL/8;
% Use an artificial value of TOL to force the program to refine.
TOL = EPS100*abs(max(sum(Q,2)));
ERR_OK = 0;
ADJUST = 1;
% Initialize storage lists before first call to SaveRectInfo.
XREFLIST = [];
QSUBLIST = zeros(0,outcls);
ADJERRLIST = zeros(0,outcls);
RECTLIST = zeros(5,0,outcls);
NRECTS = 0;
MINRECTWARN = false;
MAXNFEWARN = false;
SaveRectInfo(Qsub,esub,thetaL,thetaR,phiB,phiT);
while NRECTS > 0 && ERRBND > TOL
    % Get entries from RECTLIST corresponding to the biggest (adjusted)
    % error.
    [q,e,thetaL,thetaR,phiB,phiT,adjerr] = NextEntry;
    % Approximate integral over four subrectangles.
    [Qsub,esub] = tensor(thetaL,thetaR,phiB,phiT,q,adjerr);
    % Saved in RECTLIST is "e", a conservative estimate of the error in the
    % approximation "q" of the integral over a rectangle. Newq = sum(Qsub)
    % is a much better approximation to the integral. It is used here to
    % estimate the error in "q" and thereby determine that the estimator is
    % conservative by a factor of "ADJUST". This factor is applied to the
    % estimate of the error in "Newq" to get a more realistic estimate.
    % This scheme loses the sign of the error, so a conservative local test
    % is used to decide convergence.
    if size(Qsub,2) == 1
        % The rectangle was not subdivided because doing so would have
        % forced an evaluation on the boundary of the region of
        % integration.
        MINRECTWARN = true;
        % Move the contribution of adjerr to ERR_OK. It has already been
        % removed from ADJERRLIST.
        ERR_OK = ERR_OK + adjerr;
        ERRBND = ERR_OK + sum(ADJERRLIST);
    else
        Newq = sum(Qsub,2);
        ADJUST = min(1,max(abs(q - Newq)./e));
        Q = Q + (Newq - q);
        TOL = max(atold8,rtold8*max(abs(sum(Q,2))));
        SaveRectInfo(Qsub,esub,thetaL,thetaR,phiB,phiT);
        if NFE >= MaxFunEvals
            MAXNFEWARN = true;
            break
        end
    end
end % while
Q = sum(Q,2);

%==Nested functions======================================================

    function [Qsub,esub] = tensor(thetaL,thetaR,phiB,phiT,Qsub,esub)
        % Compute the integral with respect to theta from thetaL to thetaR
        % of the integral with respect to phi from phiB to phiT of F in
        % four blocks.
        % On the first call:
        %     The Qsub and esub input arguments are ignored. The input
        %     rectangle is subdivided into 4 subrectangles. The Qsub and
        %     esub output arguments are 4-element vectors of integral and
        %     error estimates, respectively. The first call also carries
        %     out some one-time error checking on the vectorization of F,
        %     phiBvar, and phiTvar functions.
        % On subsequent calls:
        %     Input Qsub should be the integral estimate for the input
        %     rectangle from RECTLIST, and esub should be the adjusted
        %     error estimate for this rectangle from ADJERRLIST. Normally
        %     the rectangle is subdivided, and Qsub and esub will be
        %     4-element vectors of integral and error estimates,
        %     respectively. However, if, due to roundoff error, subdividing
        %     the interval would force an evaluation on the boundary of the
        %     region, the rectangle is not subdivided, and scalars Qsub and
        %     esub are returned unchanged.
        dtheta = thetaR - thetaL;
        theta = thetaL + NARRAY*dtheta;
        if SINGULAR
            x = 0.5*(B + A) + 0.5*(B - A)*cos(theta);
            if ~FIRSTFUNEVAL && (x(1) == B || x(end) == A)
                return
            end
        else
            x = theta;
            if ~FIRSTFUNEVAL && (x(1) == A || x(end) == B)
                return
            end
        end
        X = x(ONEVEC,:);
        bottom = phiBvar(x);
        top = phiTvar(x);
        dydt = top - bottom;
        dphi = phiT - phiB;
        phi = phiB + NARRAY(:)*dphi;
        if SINGULAR
            Y = bottom(ONEVEC,:) + (0.5 + 0.5*cos(phi))*dydt;
            if ~FIRSTFUNEVAL && any(Y(1,:) == top | Y(end,:) == bottom)
                return
            end
        else
            Y = bottom(ONEVEC,:) + phi*dydt;
            if ~FIRSTFUNEVAL && any(Y(1,:) == bottom | Y(end,:) == top)
                return
            end
        end
        Z = FUN(X,Y);  NFE = NFE + 1;
        % We assume that FUN returns pages of output, so size(Z) =
        % [size(X), NUMOUTPUTS].
        if FIRSTFUNEVAL
            % Check that FUN is properly vectorized. This is important here
            % because we (otherwise) always pass in square matrices, which
            % reduces the probability of the user generating an error by
            % using matrix functions instead of elementwise functions.
            NUMOUTPUTS = size(Z,3);
            FIRSTFUNEVAL = false; % First evaluation only.
        end
        if SINGULAR
            % Full matrix formed as outer product:
            temp = 0.25*(B - A)*sin(phi)*(dydt .* sin(theta));
        else
            temp = dydt(ONEVEC,:);
        end
        Z = Z .* temp;
        for outputNum = 1 : NUMOUTPUTS
            Ztemp = Z(:,:,outputNum);
            Ztemp = [Ztemp(1:NNODES,:),Ztemp(NNODES+1:end,:)];
            r = (dtheta/4)*(dphi/4);
            % Kronrod 7 point formula tensor product.
            Qsub(outputNum,1:4) = (WT7 * reshape(WT7*Ztemp,NNODES,4))*r;
            % Gauss 3 point formula tensor product and difference with Qsub.
            esub(outputNum,1:4) = abs((WT3*reshape(WT3*Ztemp,NNODES,4))*r - Qsub(outputNum,:));
        end
    end % tensor

%--------------------------------------------------------------------------

    function SaveRectInfo(Qsub,esub,thetaL,thetaR,phiB,phiT)
        % Save information about subrectangles for which the integral is
        % not sufficiently accurate. The information is stored in four
        % arrays:
        %   RECTLIST   The columns of the RECTLIST matrix are [e;L;R;B;T],
        %              corresponding to the last 5 inputs of this function.
        %   QSUBLIST   List of Qsub values kept in the same order as
        %              RECTLIST.  This array may be complex.
        % ADJERRLIST   A list of adjusted errors kept in ascending order.
        %   XREFLIST   a cross-reference list: ADJERRLIST(idx) corresponds
        %              to RECTLIST(XREFLIST(idx)).
        %
        % NRECTS is the number of active entries in each of these lists.
        % This may be less than what is allocated. Unused entries are at
        % the end. Unused entries of ADJERRLIST and XREFLIST are zero.
        % Unused rows of RECTLIST may contain old data.
        dthetad2 = (thetaR - thetaL)/2;
        thetaM = thetaL + dthetad2;
        dphid2 = (phiT - phiB)/2;
        phiM = phiB + dphid2;
        localtol = TOL*dthetad2*dphid2/AREA;
        localtol = max(abs(localtol),EPS100*max(abs(sum(Qsub,2))));
        adjerr = max(ADJUST*esub,[],1);
        % Process each subrectangle, either adding it to the lists for
        % further subdivision or adding its adjusted error to ERR_OK.
        % Process subrectangle 1.
        if adjerr(1) > localtol
            AddToLists(Qsub(:,1),esub(:,1),thetaL,thetaM,phiB,phiM,adjerr(1));
        else
            ERR_OK = ERR_OK + adjerr(1);
        end
        % Process subrectangle 2.
        if adjerr(2) > localtol
            AddToLists(Qsub(:,2),esub(:,2),thetaM,thetaR,phiB,phiM,adjerr(2));
        else
            ERR_OK = ERR_OK + adjerr(2);
        end
        % Process subrectangle 3.
        if adjerr(3) > localtol
            AddToLists(Qsub(:,3),esub(:,3),thetaL,thetaM,phiM,phiT,adjerr(3));
        else
            ERR_OK = ERR_OK + adjerr(3);
        end
        % Process subrectangle 4.
        if adjerr(4) > localtol
            AddToLists(Qsub(:,4),esub(:,4),thetaM,thetaR,phiM,phiT,adjerr(4));
        else
            ERR_OK = ERR_OK + adjerr(4);
        end
        % Compute updated ERRBND.
        ERRBND = ERR_OK + sum(ADJERRLIST);
    end % SaveRectInfo

%--------------------------------------------------------------------------

    function AddToLists(q,e,L,R,B,T,adjerr)
        % Add [e,L,R,B,T] to RECTLIST, q to QSUBLIST, adjerr to the
        % ascending array ADJERRLIST, and cross reference index to
        % XREFLIST.
        if NRECTS >= numel(XREFLIST)
            growby = 64;
            XREFLIST(NRECTS+growby) = 0;
            RECTLIST(4+NUMOUTPUTS,NRECTS+growby) = 0;
            ADJERRLIST(NRECTS+growby) = 0;
            QSUBLIST(NUMOUTPUTS,NRECTS+growby) = 0;
        end
        NRECTS = NRECTS + 1;
        idx = find(adjerr<ADJERRLIST,1);
        if isempty(idx)
            idx = NRECTS;
        end
        % Insert sorted adjerr ascending into ADJERRLIST.
        ADJERRLIST(idx+1:NRECTS) = ADJERRLIST(idx:NRECTS-1);
        ADJERRLIST(idx) = adjerr;
        % Insert the cross-reference index into XREFLIST.
        XREFLIST(idx+1:NRECTS) = XREFLIST(idx:NRECTS-1);
        XREFLIST(idx) = NRECTS;
        % Save the data in RECTLIST.
        RECTLIST(:,NRECTS) = [e(:);L;R;B;T];
        QSUBLIST(:,NRECTS) = q;
    end % AddToLists

%--------------------------------------------------------------------------

    function [q,e,L,R,B,T,adjerr] = NextEntry
        % Return the next entry of RECTLIST and associated information.
        % This is normally the rectangle with the largest adjusted error,
        % but if the number of rectangles is greater than 2000 it will be
        % the one with the smallest adjusted error. This strategy tends to
        % arrest growth of the RECTLIST array.
        smallestFirst = NRECTS > 2000;
        if smallestFirst
            idx = XREFLIST(1);
            adjerr = ADJERRLIST(1);
        else
            idx = XREFLIST(NRECTS);
            adjerr = ADJERRLIST(NRECTS);
        end
        temp = RECTLIST(:,idx);
        e = temp(1:NUMOUTPUTS); L = temp(NUMOUTPUTS+1); R = temp(NUMOUTPUTS+2); B = temp(NUMOUTPUTS+3); T = temp(NUMOUTPUTS+4);
        q = QSUBLIST(:,idx);
        if idx ~= NRECTS
            % If idx doesn't correspond to the last active row of RECTLIST,
            % the idx row is overwritten by the last active row.
            RECTLIST(:,idx) = RECTLIST(:,NRECTS);
            QSUBLIST(:,idx) = QSUBLIST(:,NRECTS);
            XREFLIST(find(XREFLIST==NRECTS,1)) = idx;
        end
        if smallestFirst
            % We removed the first element of the sorted lists, so we must
            % shift the others.
            XREFLIST(1:NRECTS-1) = XREFLIST(2:NRECTS);
            ADJERRLIST(1:NRECTS-1) = ADJERRLIST(2:NRECTS);
        end
        XREFLIST(NRECTS) = 0;
        ADJERRLIST(NRECTS) = 0;
        NRECTS = NRECTS - 1;
    end

%--------------------------------------------------------------------------

    function p = validateAbsTol(x)
        p = true;
    end

%--------------------------------------------------------------------------

    function p = validateRelTol(x)
        p = true;
        USER_SUPPLIED_RELTOL = true;
    end

%--------------------------------------------------------------------------

    function p = validateSingular(x)
        p = true;
    end

%--------------------------------------------------------------------------

    function p = validateFailurePlot(x)
        p = true;
    end

%--------------------------------------------------------------------------

    function p = validateMaxFunEvals(x)
        p = true;
    end

%==========================================================================

end % quad2d
