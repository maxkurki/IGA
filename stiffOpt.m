function [K, F] = stiffOpt(U,V,INC,IEN,LM,gp,gw,Bz,D,R,dR_dx,J)
%Calculate stiffness matrix for optimization routines
nel = length(IEN(:,1));
Neq = max(max(LM));
K = zeros(Neq,Neq);
F = zeros(Neq,1);
mod = length(gp)^2;

for e = 1:nel

ni = INC(1,IEN(e,1)); %NURBS coordinates for elements
nj = INC(2,IEN(e,1));

if (U(ni+1) == U(ni) || V(nj+1) == V(nj)) %Skip to next element if element is zero size
   continue; 
end

neq = sum(LM(:,e) ~= 0);

Ke = zeros(neq,neq);
A = IEN(e,:);
Bez = Bz(A);
t = 0;

pos1 = 1;
for i = 1:length(gp) %Gauss quadrature
    for j = 1:length(gp)
        pos2 = (e-1)*mod+pos1;
        Jmod = J(pos2)*gw(i)*gw(j); %Jacobian modified by Gauss point weights
        t = t + Bez'*R(:,pos2)*gw(i)*gw(j)/4; %Add contribution to element thickness
        Ke = Ke + K_loc(IEN,dR_dx(:,2*pos2-1:2*pos2),D,e,LM)*Jmod; %Add contribtuion to element stiffness matrix
        pos1 = pos1 + 1;
    end
end

indx = nonzeros(LM(:,e));
K(indx,indx) = K(indx,indx) + Ke*t^3; %Add contribtuion to global stiffness matrix

end
end