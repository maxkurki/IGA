%Displacement of a cantilever

E = 2.1*10^11; % Young's modulus in N/m2
v = 0.3; % Poission's ratio

gp = [-sqrt(3/5),0,sqrt(3/5)]; %Gauss quadrature points
gw = [5/9,8/9,5/9]; %Guass quadrature point weights

D = hooke(1,E,v); % Constitutive matrix for plane stress

p = 2; % Base function degree
q = 2;
U = [zeros(1,p+1),1,2,3,4,ones(1,p+1)*5]; % Open knot vector
V = [zeros(1,q+1),0.5,ones(1,q+1)*1]; 
n = length(U)-(p+1); % Number of base functions
m = length(V)-(q+1); 

[INC,IEN] = connect(n,m,p,q); %Connectivity matrices

Ai = [];

for i = 1:n*m %Boundary conditions
   if INC(1,i) == 1
       Ai = [Ai;i,1;i,2];
   end
end


ID = destination(Ai,INC); %Coneectivity matrices
LM = location(ID,IEN);

%%
figure(3) %Plot initial structure
clf('reset')
hold on
Bx = kron(linspace(0,5,n)',ones(1,m)); %Uniform control point grid as nxm matrix
By = kron(ones(n,1),linspace(0,1,m));

B = zeros(n*m,3);
for i = 1:m %Control points as vector and control point weights
    for j = 1:n
        B(j+(i-1)*n,1) = Bx(j,i);
        B(j+(i-1)*n,2) = By(j,i);
        B(j+(i-1)*n,3) = 1;
    end
end


nelm = (n-2)*(m-2);
steps = 10;
X = zeros(steps,steps);
Y = zeros(steps,steps);
Z = -1*ones(steps,steps);
pos = linspace(-1,0.9999,steps);
odd = true;
for e = 1:nelm
    A = IEN(e,:);
    Bex = B(A,1);
    Bey = B(A,2);
    for j = 1:steps
        for k = 1:steps %Coordinates for B-spline curve
            
            [R, dR_dx, J] = Shape_function(pos(j), pos(k), e, p, q, B(:,1),B(:,2),B(:,3), U, V, INC, IEN);
            X(j,k) = Bex'*R;
            Y(j,k) = Bey'*R;
        end
    end
    if odd
        surf(X,Y,Z,'EdgeColor','none','FaceColor','b')
        odd = false;
    else
        surf(X,Y,Z,'EdgeColor','none','FaceColor','r')
        odd = true;
    end
    
end

plot(B(:,1),B(:,2),'ro','MarkerFaceColor','r','MarkerEdgeColor', 'k');

xlim([0,6]);
ylim([-0.5,1.5]);

%%

[INC,IEN] = connect(n,m,p,q);

Ai = [];

for i = 1:n*m
   if INC(1,i) == 1
       Ai = [Ai;i,1;i,2];
   end
end


ID = destination(Ai,INC);
LM = location(ID,IEN);
%%
[K,F] = stiff(U,V,INC,IEN,LM,gp,gw,p,q,B(:,1),B(:,2),B(:,3),D);
F(ID(2,n)) = F(ID(2,n)) - 1000;

d=solve(K,F);

%%
figure(11)
clf('reset')

hold on
nelm = (n-2)*(m-2);
steps = 10;
X = zeros(steps,steps);
Y = zeros(steps,steps);
Z = zeros(steps,steps);
pos = linspace(-1,0.9999,steps);
for e = 1:nelm
    A = IEN(e,:);
    Bex = B(A,1);
    Bey = B(A,2);
for j = 1:steps
for k = 1:steps %Coordinates for B-spline curve
    
    [R, dR_dx, J] = Shape_function(pos(j), pos(k), e, p, q, B(:,1),B(:,2),B(:,3), U, V, INC, IEN);
    X(j,k) = Bex'*R;
    Y(j,k) = Bey'*R;
    ee = zeros(3,1);
        P = ID(1,A);
    ee(1) = d(P + (P==0))'*((P~=0)'.*dR_dx(:,1));
    ee(3) = d(P + (P==0))'*((P~=0)'.*dR_dx(:,2));
        P = ID(2,A);
    ee(2) = d(P + (P==0))'*((P~=0)'.*dR_dx(:,2));
    ee(3) = ee(3) + d(P + (P==0))'*((P~=0)'.*dR_dx(:,1));
    ee(3) = ee(3)/2;
    st = D*ee;
    Z(j,k) =  sqrt(st(1)^2+st(2)^2 - st(1)*st(2) + 3*st(3)^2);
end
end
surf(X,Y,Z,'EdgeColor','none')
end
xlim([0,5.5]);
ylim([-0.3,1.3]);

%%
figure(12)
clf('reset')

hold on
nelm = (n-2)*(m-2);
steps = 10;
X = zeros(steps,steps);
Y = zeros(steps,steps);
Z = zeros(steps,steps);
pos = linspace(-1,0.9999,steps);
for e = 1:nelm
    A = IEN(e,:);
    BB = INC(:,A);
    Bex = B(A,1);
    Bey = B(A,2);
for j = 1:steps
    for k = 1:steps %Coordinates for B-spline curve
        
        [R, dR_dx, J] = Shape_function(pos(j), pos(k), e, p, q, B(:,1),B(:,2),B(:,3), U, V, INC, IEN);
        X(j,k) = Bex'*R;
        Y(j,k) = Bey'*R;
        ee = zeros(3,1);
        P = ID(1,A);
        ee(1) = d(P + (P==0))'*((P~=0)'.*dR_dx(:,1));
        ee(3) = d(P + (P==0))'*((P~=0)'.*dR_dx(:,2));
        X(j,k) = X(j,k) + d(P + (P==0))'*((P~=0)'.*R)*100000;
        P = ID(2,A);
        Y(j,k) = Y(j,k) + d(P + (P==0))'*((P~=0)'.*R)*100000;
        ee(2) = d(P + (P==0))'*((P~=0)'.*dR_dx(:,2));
        ee(3) = ee(3) + d(P + (P==0))'*((P~=0)'.*dR_dx(:,1));
        ee(3) = ee(3)/2;
        st = D*ee;
        Z(j,k) =  sqrt(st(1)^2+st(2)^2 - st(1)*st(2) + 3*st(3)^2);
    end
end
surf(X,Y,Z,'EdgeColor','none')
end
xlim([0,5.5]);
ylim([-0.3,1.3]);