function [z1p] = Z1p(a,alpha,u)
%MA_COEFFICENTS drag function the 2 drag coefficents in the MA matrix
%(only takes in single values)
%   [z1p,z2p,dz1p,dz2p] = MA_coefficents(a,alpha,u)

z1 = Z1(a,alpha);

S=size(alpha);
if S==1
    S=size(u);
end
z1p=zeros(S(1),S(2));

for ii=1:S(1)
    for jj=1:S(2)
        if size(alpha)~=1
        alpha2=alpha(ii,jj);
        a2=a(ii,jj);
        else
           alpha2=alpha;
           a2=a;
        end
        u2=u(ii,jj);
        
        k1=sqrt(1-u2);
        k2=sqrt(1+u2);
        
        if abs(alpha2 - 1) < 1e-6
            if abs(u2) < 1e-6
                zt1=7/(3*sqrt(2));
            elseif abs(u2 - 1) < 1e-6 || abs(u2 + 1) < 1e-6
                zt1=8/3;
            else
                zt1= (8*(k2-k1)+u2.*(-4*(k1+k2)+u2.*(19*(k1-k2)+3*u2.*(u2.*( k2-k1)+4*(k1+k2))))) ./ (3*(u2.^3).*sqrt(2-2*u2.^2));
            end
        else
            
            l=log(-(((1+u2-alpha2^2+sqrt((1+u2).*(1+u2-2*u2*alpha2^2+(u2-1).*alpha2^4))).*(u2-1+alpha2^2-sqrt((u2-1).*(alpha2^4-1+u2.*(alpha2^2-1)^2))))./((alpha2^2).*(2-2*alpha2^2 + (u2.^2).*(alpha2^2-2)))));
            
            if abs(u2 - 1) < 1e-6 || abs(u2 + 1) < 1e-6
                zt1=2/(-1+alpha2^2)-(2*(-2+alpha2^2)*log((1+sqrt(1-alpha2^2))/alpha2))/(1-alpha2^2)^(3/2);
            else
                zt11=-((l.*(-2+alpha2^2))/(1-alpha2^2)^(3/2));
                zt12=-(2*k1*(-(k1^2)*(k2^4)+(1+2*u2-u2^3)*alpha2^2)*sqrt(1-alpha2^4-u2*((-1+alpha2^2)^2))+2*k2*(-(k1^4)*(k2^2)+(1-2*u2+u2^3)*alpha2^2)*sqrt(1-alpha2^4+u2*((-1+alpha2^2)^2)))/((-1+alpha2^2)*(2-2*alpha2^2+(u2^2)*(-2+alpha2^2))*sqrt((1-alpha2^2)*((u2^4)*((-1+alpha2^2)^2)+((1+alpha2^2)^2)-2*(u2^2)*(1+alpha2^4))));
                zt1=zt11+zt12;
            end
        end
        
        z1p(ii,jj) = real(z1(ii,jj)-2*pi*zt1/a2);
    end
end

end

