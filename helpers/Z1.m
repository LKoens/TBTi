function [z1] = Z1(a,alpha)

    z1 = 4*pi*((alpha.^2-2).*acos(1./alpha)+(alpha.^2-1).^0.5) ./ (a .*(alpha.^2-1).^(1.5));
    % Any alpha that are too close to unity should have the output replaced
    % with the limiting value.
    mask = abs(alpha - 1) < 1e-6;
    z1(mask) = 16*pi./(3*a(mask));

end

