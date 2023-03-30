function out = B_integrand(t, s, r1, r2, r3, rho, epps, d)
%% B_INTEGRAND will return the combined integrand needed for evaluating B. The
%ordering of components matches that of unpack_B_integrals.

    out = zeros(23,1);

    % Evaluate functions once.
    r1s = r1(s);
    r2s = r2(s);
    r3s = r3(s);

    r1t = r1(t);
    r2t = r2(t);
    r3t = r3(t);

    rhos = rho(s);
    rhot = rho(t);

    denom1 = ((r1s-r1t).^2+(r2s-r2t).^2+(r3s-r3t).^2 +(epps^2) .* (rhos.^2 +rhot.^2)).^(3/2);
    denom2 = ((r1s-r1t).^2+(r2s-r2t).^2+(r3s+r3t-2 * d).^2 +(epps^2) .* (rhos.^2 +rhot.^2)).^(3/2);
    denom3 = ((r1s-r1t).^2+(r2s-r2t).^2+(r3s+r3t-2 * d).^2 +(epps^2) .* (rhos.^2 +rhot.^2)).^(5/2);

    % Stokeslet components.
    out(1) = 1 ./ sqrt((r1s-r1t).^2+(r2s-r2t).^2+(r3s-r3t).^2 + (epps^2) .* (rhos.^2 +rhot.^2));
    out(2) = (r1s-r1t).^2 ./ denom1;
    out(3) = (r1s-r1t) .* (r2s-r2t) ./ denom1;
    out(4) = (r1s-r1t) .* (r3s-r3t) ./ denom1;
    out(5) = (r2s-r2t) .* (r2s-r2t) ./ denom1;
    out(6) = (r2s-r2t) .* (r3s-r3t) ./ denom1;
    out(7) = (r3s-r3t) .* (r3s-r3t) ./ denom1;

    % Stokeslet boundary corrections.
    out(8)  = 1 ./ sqrt((r1s-r1t).^2+(r2s-r2t).^2+(r3s+r3t-2*d).^2 + (epps^2) .* (rhos.^2 +rhot.^2));
    out(9)  = (r1s-r1t) .* (r1s-r1t) ./ denom2;
    out(10) = (r1s-r1t) .* (r2s-r2t) ./ denom2;
    out(11) = (r1s-r1t) .* (r3s+r3t-2 * d) ./ denom2;
    out(12) = (r2s-r2t) .* (r2s-r2t) ./ denom2;
    out(13) = (r2s-r2t) .* (r3s+r3t-2 * d) ./ denom2;
    out(14) = (r3s+r3t-2 * d) .* (r3s+r3t-2 * d) ./ denom2;

    % Source dipole. %THE H0 STUFF HERE LOOKS TO BE WRONG. IN BENS FORMULA
    % BOTH ARE INTEGRATED OVER.
    out(15) = (r3t-d) .* (r3s-d) ./ denom2;
    out(16) = (r3t-d) .* (r3s-d) .* (r1s-r1t) .* (r1s-r1t) ./ denom3;
    out(17) = (r3t-d) .* (r3s-d) .* (r1s-r1t) .* (r2s-r2t) ./ denom3;
    out(18) = (r3t-d) .* (r3s-d) .* (r1s-r1t) .* (r3s+r3t-2*d) ./ denom3;
    out(19) = (r3t-d) .* (r3s-d) .* (r2s-r2t) .* (r2s-r2t) ./ denom3;
    out(20) = (r3t-d) .* (r3s-d) .* (r2s-r2t) .* (r3s+r3t-2*d) ./ denom3;
    out(21) = (r3t-d) .* (r3s-d) .* (r3s+r3t-2 * d) .* (r3s+r3t-2*d) ./ denom3;

    % Stresslet. %THE H0 STUFF HERE LOOKS TO BE WRONG
    out(22) = (r3t-d) .* (r1s-r1t) ./ denom2;
    out(23) = (r3t-d) .* (r2s-r2t) ./ denom2;

end