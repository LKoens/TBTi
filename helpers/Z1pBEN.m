function [z1p] = Z1pBEN(a,alpha,u)
%MA_COEFFICENTS drag function the 2 drag coefficents in the MA matrix
%(only takes in single values)
%   [z1p,z2p,dz1p,dz2p] = MA_coefficents(a,alpha,u)

beta = alpha.^2 - 1;
rootBeta = beta.^0.5;

Z1prime = 4*pi*((beta - 1).*acos(1./alpha) + rootBeta) ./ (a .* rootBeta.^3);


        
        k1=sqrt(1-u2);
        k2=sqrt(1+u2);
        
        if alpha2==1%abs(alpha2-1) < 1e-9
            z1= 8/(3*a2);
            if u2==0
                zt1=7/(3*sqrt(2));
            elseif u2^2==1
                zt1=8/3;
            else
                zt1= (8*(k2-k1)+u2.*(-4*(k1+k2)+u2.*(19*(k1-k2)+3*u2.*(u2.*( k2-k1)+4*(k1+k2))))) ./ (3*(u2.^3).*sqrt(2-2*u2.^2));
            end
        else
            z1=2*((alpha2^2-2).*acos(1/alpha2)+sqrt(alpha2^2-1))/(a2*(alpha2^2-1)^(3/2));
            
            l=log(-(((1+u2-alpha2^2+sqrt((1+u2).*(1+u2-2*u2*alpha2^2+(u2-1).*alpha2^4))).*(u2-1+alpha2^2-sqrt((u2-1).*(alpha2^4-1+u2.*(alpha2^2-1)^2))))./((alpha2^2).*(2-2*alpha2^2 + (u2.^2).*(alpha2^2-2)))));
            
            if u2^2==1
                zt1=2/(-1+alpha2^2)-(2*(-2+alpha2^2)*log((1+sqrt(1-alpha2^2))/alpha2))/(1-alpha2^2)^(3/2);
            else
                zt11=-((l.*(-2+alpha2^2))/(1-alpha2^2)^(3/2));
                zt12=-(2*k1*(-(k1^2)*(k2^4)+(1+2*u2-u2^3)*alpha2^2)*sqrt(1-alpha2^4-u2*((-1+alpha2^2)^2))+2*k2*(-(k1^4)*(k2^2)+(1-2*u2+u2^3)*alpha2^2)*sqrt(1-alpha2^4+u2*((-1+alpha2^2)^2)))/((-1+alpha2^2)*(2-2*alpha2^2+(u2^2)*(-2+alpha2^2))*sqrt((1-alpha2^2)*((u2^4)*((-1+alpha2^2)^2)+((1+alpha2^2)^2)-2*(u2^2)*(1+alpha2^4))));
                zt1=zt11+zt12;
            end
        end
        
        z1p(ii,jj) =real(2*pi*(z1-zt1/a2));
    end
end

end

