function out = c4Integrand(s,phi,st,pt,u,sFun,S1,S2,S3,S1e,S2e,S3e)
    out = zeros([size(s),7]);

    % Precompute commonly used values.
    S1t = S1(st,pt);
    S2t = S2(st,pt);
    S3t = S3(st,pt);

    S1v = S1(s,phi);
    S2v = S2(s,phi);
    S3v = S3(s,phi);

    S1eu = S1e(u,pt);
    S2eu = S2e(u,pt);
    S3eu = S3e(u,pt);

    S1e3 = S1e(sFun(s),phi);
    S2e3 = S2e(sFun(s),phi);
    S3e3 = S3e(sFun(s),phi);

    denom1 = ((S1t-S1v).^2+(S2t-S2v).^2+(S3t-S3v).^2).^0.5;
    denom3 = denom1.^3;

    out(:,:,1) = 1 ./ denom1;
    out(:,:,2) = (S1t-S1v) .* (S1t-S1v) ./ denom3;
    out(:,:,3) = (S1t-S1v) .* (S2t-S2v) ./ denom3;
    out(:,:,4) = (S1t-S1v) .* (S3t-S3v) ./ denom3;
    out(:,:,5) = (S2t-S2v) .* (S2t-S2v) ./ denom3;
    out(:,:,6) = (S2t-S2v) .* (S3t-S3v) ./ denom3;
    out(:,:,7) = (S3t-S3v) .* (S3t-S3v) ./ denom3;
end