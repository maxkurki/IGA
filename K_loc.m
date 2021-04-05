function [Ke] = K_loc(IEN,dR_dx,D,e,LM)
%Calculate element stiffness matrix

nen = length(IEN(1,:));

B = zeros(3,2*nen);

for i = 1:nen %B-matrix from basis function derivatives
   B(:,2*i-1:2*i) = [dR_dx(i,1),0;0,dR_dx(i,2);dR_dx(i,2),dR_dx(i,1)]; 
end

Ke = B'*D*B; %Element stiffness matrix


Ke = Ke + 1; %Add one to avoid zeros in stiffness matrix

for p = 1:nen*2 %Set entries corresponding to fixed control points to zero 
   if (LM(p,e) == 0)
      Ke(:,p) = 0;
      Ke(p,:) = 0;
   end
end

Ke = nonzeros(Ke); %Reshape stiffness matrix without entries corresponding to fixed control points
Ke = Ke - 1;
neq = sum(LM(:,e) ~= 0);

Ke = reshape(Ke,neq,neq);
end



