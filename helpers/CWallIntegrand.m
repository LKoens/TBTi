function out = CWallIntegrand(s,phi,st,pt,S1,S2,S3,S1e,S2e,S3e)
%% CWALLINTEGRAND will return the combined integrand needed for evaluating
%the contribution of the wall images.

    out = zeros([size(s),16]);

    % Precompute commonly used values.
    S1t = S1(st,pt);
    S2t = S2(st,pt);
    S3t = S3(st,pt);

    S1v = S1(s,phi);
    S2v = S2(s,phi);
    S3v = S3(s,phi);

    denom1 = ((S1t-S1v).^2+(S2t-S2v).^2+(S3t+S3v).^2).^0.5;
    denom3 = denom1.^3;
    denom5 = denom1.^5;

    % Mirror stokeslet components
    out(:,:,1) =  1 ./ denom1;
    out(:,:,2) = (S1t-S1v) .* (S1t-S1v) ./ denom3;
    out(:,:,3) = (S1t-S1v) .* (S2t-S2v) ./ denom3;
    out(:,:,4) = (S1t-S1v) .* (S3t+S3v) ./ denom3;
    out(:,:,5) = (S2t-S2v) .* (S2t-S2v) ./ denom3;
    out(:,:,6) = (S2t-S2v) .* (S3t+S3v) ./ denom3;
    out(:,:,7) = (S3t+S3v) .* (S3t+S3v) ./ denom3;
    
    % Mirror source dipole
    out(:,:,8)  = S3v .* S3t ./ denom3;
    out(:,:,9)  = S3v .* S3t .* (S1t-S1v) .* (S1t-S1v) ./ denom5;
    out(:,:,10) = S3v .* S3t .* (S1t-S1v) .* (S2t-S2v) ./ denom5;
    out(:,:,11) = S3v .* S3t .* (S1t-S1v) .* (S3t+S3v) ./ denom5;
    out(:,:,12) = S3v .* S3t .* (S2t-S2v) .* (S2t-S2v) ./ denom5;
    out(:,:,13) = S3v .* S3t .* (S2t-S2v) .* (S3t+S3v) ./ denom5;
    out(:,:,14) = S3v .* S3t .* (S3t+S3v) .* (S3t+S3v) ./ denom5;

    % Mirror Stresslet
    out(:,:,15) = S3v .* (S1t-S1v) ./ denom3;
    out(:,:,16) = S3v .* (S2t-S2v) ./ denom3;

end