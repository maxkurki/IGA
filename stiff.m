function [K, F] = stiff(U,V,INC,IEN,LM,gp,gw,p,q,Bx,By,w,D)
%Calculate stiffness matrix, same as stiffOpt but without precalculated
%basis functions values, derivates, and jacobians

nel = length(IEN(:,1));
Neq = max(max(LM));
K = zeros(Neq,Neq);
F = zeros(Neq,1);

for e = 1:nel

ni = INC(1,IEN(e,1)); %NURBS coordinates for elements
nj = INC(2,IEN(e,1));

if (U(ni+1) == U(ni) || V(nj+1) == V(nj))
   continue; 
end

neq = sum(LM(:,e) ~= 0);

Ke = zeros(neq,neq);
Fe = zeros(neq,1);

for i = 1:length(gp)
    for j = 1:length(gp)
        [R,dR_dx,J] = Shape_function(gp(i),gp(j),e,p,q,Bx,By,w,U,V,INC,IEN);
        Jmod = J*gw(i)*gw(j);
        Ke = Ke + K_loc(IEN,dR_dx,D,e,LM)*Jmod;
    end
end

indx = nonzeros(LM(:,e));
K(indx,indx) = K(indx,indx) + Ke;

end
end

