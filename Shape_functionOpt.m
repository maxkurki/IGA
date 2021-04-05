function [R, dR_dx, J1, J2] = Shape_functionOpt(xi_t, eta_t, e, p, q, Bx,By,W, U, V, INC, IEN)
%Adapted from pp. 100-102 in Hughes
%Shape function for IGA element

nen = (p+1)*(q+1); %Number of NURBS basis functions
R = zeros(nen,1); %Array to store NURBS basis function values
dR_dx = zeros(nen,2); %Array to store NURBS basis function derivatives

dR_dxi = zeros(nen,2);   %Array to store NURBS basis funtcion derivates with 
                         %with respect to parametric coordinates
dx_dxi = zeros(2,2); %Derivative of physical coordinates with respect to 
                     %parametric coordinates

dxi_dtxi = zeros(2,2); %Derivatives of parametric coordinates with 
                       %with respect to parent element coordinates

loc_num = 0; %Local basis function number

sum_xi = 0;   %Sums for calculating NURBS functions values
sum_eta = 0;  %and derivatives
sum_tot = 0;

ni = INC(1,IEN(e,1)); %NURBS coordinates for elements
nj = INC(2,IEN(e,1));

A = IEN(e,:);

xi = (U(ni)+U(ni+1))/2 + (U(ni+1)-U(ni))*xi_t/2;  %Parametric coordinates from Gauss points
eta = (V(nj)+V(nj+1))/2 + (V(nj+1)-V(nj))*eta_t/2;%in parent element

n = length(U) - p - 1;
m = length(V) - q - 1;

N = BasisFuns2(xi,p,n,1,U);  %Calculate basis function values and derivatives 
M = BasisFuns2(eta,q,m,1,V); %at the parametric coordinates

for j = 0:q      %Numerators and denominators for NURBS basis functions
    for i = 0:p  %and derivatives
        loc_num = loc_num + 1;
        R(loc_num) = N(p+1-i,1)*M(q+1-j,1)*W(A(loc_num));
        sum_tot = sum_tot + R(loc_num);
        dR_dxi(loc_num,1) = N(p+1-i,2)*M(q+1-j,1)*W(A(loc_num));
        sum_xi = sum_xi + dR_dxi(loc_num,1);
        dR_dxi(loc_num,2) = N(p+1-i,1)*M(q+1-j,2)*W(A(loc_num));
        sum_eta = sum_eta + dR_dxi(loc_num,2);
    end
end

for i = 1:nen %Finalizing NURBS basis functions and derivatives
   R(i) = R(i)/sum_tot;
   dR_dxi(i,1) = (dR_dxi(i,1)*sum_tot - R(i)*sum_xi)/(sum_tot^2);
   dR_dxi(i,2) = (dR_dxi(i,2)*sum_tot - R(i)*sum_eta)/(sum_tot^2);
end

loc_num = 0;

for j = 0:q     %Gradient of mapping from parametric coordinates to 
    for i = 0:p %physical coordinates
        loc_num = loc_num + 1;
        
        for jj = 1:2
            dx_dxi(1,jj) = dx_dxi(1,jj) + Bx(A(loc_num))*dR_dxi(loc_num,jj);
            dx_dxi(2,jj) = dx_dxi(2,jj) + By(A(loc_num))*dR_dxi(loc_num,jj);
        end
    end
end
dxi_dx = inv(dx_dxi); %Inverse of above

for i = 1:nen   %Derivatives of NURBS basis functions with respect to
   for aa = 1:2 %physical coordinates
       for bb = 1:2
           dR_dx(i,aa) = dR_dx(i,aa) + dR_dxi(i,bb)*dxi_dx(bb,aa);
       end
   end
end

dxi_dtxi(1,1) = (U(ni+1)-U(ni))/2; %Gradient of mapping from parent element coordinates
dxi_dtxi(2,2) = (V(nj+1)-V(nj))/2; %to parametric coordinates

J1 = det(dx_dxi); %Jacobian determinant of mapping from parametric coordinates to physical coordinates
J2 = det(dxi_dtxi); %Jacobian determinant of mapping from parent element coordinates to parametric coordinates
end

