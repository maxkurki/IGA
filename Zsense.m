function [dg0,dg1,vol] = Zsense(d,U,V,INC,IEN,LM,gp,gw,Bz,D,R,dR_dx,J)
% Calculates sensitivity of objective function and volume constraint for
% topology optimization as well as the volume of the structure
nel = length(IEN(:,1));
dg0 = zeros(length(INC),1);
dg1 = zeros(length(INC),1);
vol = 0;
mod = length(gp)^2;

for e = 1:nel

ni = INC(1,IEN(e,1)); %NURBS coordinates for elements
nj = INC(2,IEN(e,1));

if (U(ni+1) == U(ni) || V(nj+1) == V(nj)) %Skip to next element if element is zero size
   continue; 
end

neq = sum(LM(:,e) ~= 0);
Ke = zeros(neq,neq);
indx = nonzeros(LM(:,e));
A = IEN(e,:); %Global shape functions numbers for current element
Bez = Bz(A); %Control points associated with current element
dtdB = zeros(length(A),1);
area = 0;
t = 0;

pos1 = 1;
for i = 1:length(gp) %Gauss quadrature
    for j = 1:length(gp)
        pos2 = (e-1)*mod+pos1;
        Jmod = J(pos2)*gw(i)*gw(j); %Jacobian modified by Gauss point weights
        t = t + Bez'*R(:,pos2)*gw(i)*gw(j)/4; %Add contribution to element thickness
        dtdB = dtdB + R(:,pos2)*gw(i)*gw(j)/4; %Add contribution to derivative of thickness with respect to control points
        area = area + Jmod; %Add contribution to element area
        Ke = Ke + K_loc(IEN,dR_dx(:,2*pos2-1:2*pos2),D,e,LM)*Jmod; %Add contribution to element stiffness matrix
        pos1 = pos1 + 1;
    end
end

vol = vol + area*t; %Add volume of current element to total volume

for k = 1:length(A)
    dg0(A(k)) = dg0(A(k)) -d(indx)'*Ke*d(indx)*3*t^2*dtdB(k); %Add contribtuion to sensitivities
    dg1(A(k)) = dg1(A(k)) + area*dtdB(k);
end

end
end
