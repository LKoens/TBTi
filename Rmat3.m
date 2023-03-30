function [R] = Rmat3(SO0,SNO,numS,numPhi)
%RMATv2 Summary of this function goes here
%   Only calculates force

R=zeros(6,6);

ds = 2/(numS); %distance between consecutive s points.
dphi = 2*pi/(numPhi);

V0 = zeros(3*(numPhi)*numS,6);
for i=1:numS
    for j=1:(numPhi)
        V0((1:3)+3*(i-1)+3*numS*(j-1),1)=8*pi*[1,0,0];
        V0((1:3)+3*(i-1)+3*numS*(j-1),2)=8*pi*[0,1,0];
        V0((1:3)+3*(i-1)+3*numS*(j-1),3)=8*pi*[0,0,1];
    end
end


f0=(eye(3*(numPhi)*numS) + SNO)\SO0*V0; %solution using geometric series

% Compute R using f_0.
R(1:3,:) = squeeze(sum(reshape(f0,3,[],6),2)) * dphi * ds; % Total force.

end

