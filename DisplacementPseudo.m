%Displacement of a two-element cantilever

E = 2.1*10^11; % Young's modulus in N/m2
v = 0.3; % Poission's ratio

gp = [-sqrt(3/5),0,sqrt(3/5)]; %Gauss quadrature points
gw = [5/9,8/9,5/9]; %Guass quadrature point weights

D = hooke(1,E,v); % Constitutive matrix for plane stress

p = 2; % Base function degree
U = [zeros(1,p+1),1,ones(1,p+1)*2]; % Open knot vector
q = 2;
V = [zeros(1,q+1),ones(1,q+1)*1];  
n = length(U)-(p+1); % Number of base functions
m = length(V)-(q+1); 

[INC,IEN] = connect(n,m,p,q); %Connectivity matrices

Ai = [];

for i = 1:n*m %Boundary conditions
   if INC(1,i) == 1
       Ai = [Ai;i,1;i,2];
   end
end


ID = destination(Ai,INC); %Connectivity matrices
LM = location(ID,IEN);



%%

Bx = kron(linspace(0,30,n)',ones(1,m)); %Uniform control point grid as nxm matrix
By = kron(ones(n,1),linspace(0,10,m));

B = zeros(n*m,3);
for i = 1:m %Control points as vector and control point weights
    for j = 1:n
        B(j+(i-1)*4,1) = Bx(j,i);
        B(j+(i-1)*4,2) = By(j,i);
        B(j+(i-1)*4,3) = 1;
    end
end

figure(3) %Plot initial structure
clf('reset')
hold on
nelm = (n-2)*(m-2); %Number of elements
steps = 10; %Evaluation points in each element
X = zeros(steps,steps);
Y = zeros(steps,steps);
Z = -1*ones(steps,steps);
pos = linspace(-0.9999,0.9999,steps); %Parent element domain is [-1,1]^2
odd = true;
for e = 1:nelm
    A = IEN(e,:); %Global shape functions numbers for current element
    Bex = B(A,1); %Control points associated with current element
    Bey = B(A,2);
    for j = 1:steps
        for k = 1:steps %Coordinates for B-spline curve
            
            [R, ~, ~] = Shape_function(pos(j), pos(k), e, p, q, B(:,1),B(:,2),B(:,3), U, V, INC, IEN); %Basis function values in current evaluation point
            X(j,k) = Bex'*R; %Coordinates resulting from current evaluation point
            Y(j,k) = Bey'*R;
        end
    end
    if odd
        surf(X,Y,Z,'EdgeColor','none','FaceColor','b') %Plot surface in alternating color for visibility
        odd = false;
    else
        surf(X,Y,Z,'EdgeColor','none','FaceColor','r')
        odd = true;
    end
    
end

plot(B(:,1),B(:,2),'ro','MarkerFaceColor','r','MarkerEdgeColor', 'k'); %Plot control points

xlim([0,33]);
ylim([-2,12]);

%%
[K,F] = stiff(U,V,INC,IEN,LM,gp,gw,p,q,B(:,1),B(:,2),B(:,3),D); %Stiffness matrix and force vector
F(ID(2,n)) = F(ID(2,n)) - 100000;

d=solve(K,F); %Solve for displacements

%%
figure(16) %Plot displacement field
clf('reset')

hold on
nelm = (n-2)*(m-2); %Number of elements
steps = 10; %Evaluation points in each element
X = zeros(steps,steps);
Y = zeros(steps,steps);
Z = -1*ones(steps,steps);
pos = linspace(-0.9999,0.9999,steps); %Parent element domain is [-1,1]^2
for e = 1:nelm
    A = IEN(e,:); %Global shape functions numbers for current element
    Bex = B(A,1); %Control points associated with current element
    Bey = B(A,2);
    for j = 1:steps
        for k = 1:steps %Coordinates for B-spline curve
            
            [R, dR_dx, ~] = Shape_function(pos(j), pos(k), e, p, q, B(:,1),B(:,2),B(:,3), U, V, INC, IEN); %Basis functions values and derivatives in current evaluation point
            X(j,k) = Bex'*R; %Coordinates resulting from current evaluation point
            Y(j,k) = Bey'*R;
            ee = zeros(3,1);
            P = ID(1,A); %Index of displacements in x-drection associated with current element
            X(j,k) = X(j,k) + d(P + (P==0))'*((P~=0)'.*R)*1000; %Displaced x-position with magnification of 1000
            ee(1) = d(P + (P==0))'*((P~=0)'.*dR_dx(:,1)); %Strain exx
            ee(3) = d(P + (P==0))'*((P~=0)'.*dR_dx(:,2)); %Strain exy
            P = ID(2,A); %Index of displacements in y-drection associated with current element
            Y(j,k) = Y(j,k) + d(P + (P==0))'*((P~=0)'.*R)*1000; %Displaced y-position with magnification of 1000
            ee(2) = d(P + (P==0))'*((P~=0)'.*dR_dx(:,2)); %Strain eyy
            ee(3) = ee(3) + d(P + (P==0))'*((P~=0)'.*dR_dx(:,1)); %Strain exy
            ee(3) = ee(3)/2; %Strain exy
            st = D*ee; %Stress
            Z(j,k) =  sqrt(st(1)^2+st(2)^2 - st(1)*st(2) + 3*st(3)^2); %Von mises stress
        end
    end
    surf(X,Y,Z,'EdgeColor','none') %Plot surface
end

xlim([0,33]);
ylim([-2,12]);