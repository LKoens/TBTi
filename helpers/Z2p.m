function [z2p] = Z2p(a,alpha,u)
%MA_COEFFICENTS drag function the 2 drag coefficents in the MA' matrix
%(only takes in single values)
%   [z1p,z2p,dz1p,dz2p] = MA_coefficents(a,alpha,u)

z2 = Z2(a,alpha);

S=size(alpha);
if S==1
    S=size(u);
end
z2p=zeros(S(1),S(2));

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
                zt2=sqrt(2);
            elseif abs(u2 - 1) < 1e-6 || abs(u2 + 1) < 1e-6
                zt2=2;
            else
                zt2=sqrt(2)*(k2-k1)./u2;
            end
        else
            
            l=log(-(((1+u2-alpha2^2+sqrt((1+u2).*(1+u2-2*u2*alpha2^2+(u2-1).*alpha2^4))).*(u2-1+alpha2^2-sqrt((u2-1).*(alpha2^4-1+u2.*(alpha2^2-1)^2))))./((alpha2^2).*(2-2*alpha2^2 + (u2.^2).*(alpha2^2-2)))));
            
            zt2=l./sqrt(1-alpha2^2);
        end
        
        z2p(ii,jj) =real(z2(ii,jj)-2*pi*zt2/a2);
        
    end
end
end

