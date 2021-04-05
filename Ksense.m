function [Ke] = Ksense(IEN,dR_dx,dG_da,J,dJ_da,D,e,LM)
%Calculate sensitivity of element stiffness matrix

nen = length(IEN(1,:));

B = zeros(3,2*nen);

for i = 1:nen %B-matrix for basis function derivates
   B(:,2*i-1:2*i) = [dR_dx(i,1),0;0,dR_dx(i,2);dR_dx(i,2),dR_dx(i,1)]; 
end

dB_da = zeros(3,2*nen); %B-matrix derivative with respect to design variable

for i = 1:nen
   dB_da(:,2*i-1:2*i) = [dG_da(i,1),0;0,dG_da(i,2);dG_da(i,2),dG_da(i,1)]; 
end

Ke = dB_da'*D*B*J + B'*D*dB_da*J + B'*D*B*dJ_da; %Element stiffness matrix derivative with respect to design variable


Ke = Ke + 1; %Same as K_loc, reshape without entries corresponding to fixed control points

for p = 1:nen*2
   if (LM(p,e) == 0)
      Ke(:,p) = 0;
      Ke(p,:) = 0;
   end
end

Ke = nonzeros(Ke);
Ke = Ke - 1;
neq = sum(LM(:,e) ~= 0);

Ke = reshape(Ke,neq,neq);
end

