E = 2.1*10^11; % Young's modulus in N/m2
v = 0.3; % Poission's ratio

gp = [-sqrt(3/5),0,sqrt(3/5)]; %Gauss quadrature points
gw = [5/9,8/9,5/9]; %Guass quadrature point weights

D = hooke(1,E,v); % Constitutive matrix for plane stress

p = 2; % Base function degree
U = [zeros(1,p),linspace(0,2,36),ones(1,p)*2]; % Open knot vector
q = 2;
V = [zeros(1,q),linspace(0,1,12),ones(1,q)*1]; 
n = length(U)-(p+1); % Number of base functions
m = length(V)-(q+1); 

[INC,IEN] = connect(n,m,p,q); %Connectivity matrices

Ai = [n,2];

for i = 1:n*m %Boundary conditions
   if INC(1,i) == 1
       Ai = [Ai;i,1];
   end
end


ID = destination(Ai,INC); %Connectivity matrices
LM = location(ID,IEN);

Bx = kron([0,linspace(0.1,29.9,n-2),30]',ones(1,m)); %Control point grid as nxm matrix
By = kron(ones(n,1),[0,linspace(0.1,9.9,m-2),10]);

B = zeros(n*m,3);
for i = 1:m %Control points as vector and control point weights
    for j = 1:n
        B(j+(i-1)*n,1) = Bx(j,i);
        B(j+(i-1)*n,2) = By(j,i);
        B(j+(i-1)*n,3) = 1;
    end
end

nelm = (n-p)*(m-q);

row = (p+1)*(q+1); 
col = nelm*length(gp)^2;
R = zeros(row,col);
dR_dx = zeros(row,2*col);
J = zeros(1,col);
mod = length(gp)^2;
for e = 1:nelm %Precalculating shape function values, derivatives, and jacobians in all Gauss points for all elements
    pos1 = 1;
    for i = 1:length(gp)
        for j = 1:length(gp)
            pos2 = (e-1)*mod+pos1;
            [R(:,pos2),dR_dx(:,2*pos2-1:2*pos2),J(pos2)] = Shape_function(gp(i),gp(j),e,p,q,B(:,1),B(:,2),B(:,3),U,V,INC,IEN);
            pos1 = pos1 + 1;
        end
    end
end

Vmax = 0.35*10*30; % Maximum allowed volume

T = 0.35*ones(length(INC),1); % Initial element thickness;
T_old1 = T; %MMA parameters
T_old2 = T;
Tmin = 0.001;
Tmax = 1;
low = Tmin;
upp = Tmax;
C_old = 0;

run = true;
iter = 0;
while run
    iter = iter + 1;
    
    [K,F] = stiffOpt(U,V,INC,IEN,LM,gp,gw,p,q,B(:,1),B(:,2),T,B(:,3),D,R,dR_dx,J); %Stiffness matrix and force vector
    F(ID(2,n*(m-1)+1)) = F(ID(2,n*(m-1)+1)) - 10000;
    
    d=solve(K,F); %Solve for displacements
    [dg0,dg1,vol] = Zsense(d,U,V,INC,IEN,LM,gp,gw,p,q,B(:,1),B(:,2),T,B(:,3),D,R,dR_dx,J); %Sensitivities and volume
    fval = vol-Vmax; %current value of constraint function
    C1 = d'*K'*d; %Compliance
        [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp] = ... %Use MMA to clalculate new values
        mmasub(1,length(INC),iter,T,Tmin,Tmax,T_old1,T_old2, ...
        C1,dg0,fval,dg1,low,upp,1,0,1000,1);
    T_old2 = T_old1; %Update values
    T_old1 = T;
    T = xmma;
    if norm(C1-C_old)<0.001 %Check convergence
        run = false;
    end
    C_old = C1;
end

%%

figure(23) %Plot optimized structure
clf('reset')
hold on
nelm = (n-p)*(m-q); %number of elements
steps = 10; %Evaluation points in each element
X = zeros(steps,steps);
Y = zeros(steps,steps);
Z = zeros(steps,steps);
pos = linspace(-0.9999,0.9999,steps); %Parent element domain is [-1,1]^2
for e = 1:nelm
    A = IEN(e,:); %Global shape functions numbers for current element
    Bex = B(A,1); %Control points associated with current element
    Bey = B(A,2);
    Bez = T(A);
for j = 1:steps
    for k = 1:steps %Coordinates for B-spline curve  
        [R, ~, ~] = Shape_function(pos(j), pos(k), e, p, q, B(:,1),B(:,2),B(:,3), U, V, INC, IEN); %Basis function values in current evaluation point
        X(j,k) = Bex'*R; %Coordinates resulting from current evaluation point
        Y(j,k) = Bey'*R;
        Z(j,k) = Bez'*R;
        
    end
end
surf(X,Y,Z,'EdgeColor','none') %Plot surface
end
xlim([0,33]);
ylim([-2,12]);

%%
%Same as above but with p=q=4
p = 4; %
U = [zeros(1,p),linspace(0,2,36),ones(1,p)*2]; 
q = 4;
V = [zeros(1,q),linspace(0,1,12),ones(1,q)*1]; 
n = length(U)-(p+1); 
m = length(V)-(q+1); 

[INC,IEN] = connect(n,m,p,q);

Ai = [n,2];

for i = 1:n*m
   if INC(1,i) == 1
       Ai = [Ai;i,1];
   end
end


ID = destination(Ai,INC);
LM = location(ID,IEN);

Bx = kron([0,linspace(0.1,29.9,n-2),30]',ones(1,m));
By = kron(ones(n,1),[0,linspace(0.1,9.9,m-2),10]);

B = zeros(n*m,3);
for i = 1:m
    for j = 1:n
        B(j+(i-1)*n,1) = Bx(j,i);
        B(j+(i-1)*n,2) = By(j,i);
        B(j+(i-1)*n,3) = 1;
    end
end

nelm = (n-p)*(m-q);

row = (p+1)*(q+1);
col = nelm*length(gp)^2;
R = zeros(row,col);
dR_dx = zeros(row,2*col);
J = zeros(1,col);
mod = length(gp)^2;
for e = 1:nelm
    pos1 = 1;
    for i = 1:length(gp)
        for j = 1:length(gp)
            pos2 = (e-1)*mod+pos1;
            [R(:,pos2),dR_dx(:,2*pos2-1:2*pos2),J(pos2)] = Shape_function(gp(i),gp(j),e,p,q,B(:,1),B(:,2),B(:,3),U,V,INC,IEN);
            pos1 = pos1 + 1;
        end
    end
end

run = true;

Vmax = 0.35*10*30; 


T = 0.35*ones(length(INC),1);
T_old1 = T; 
T_old2 = T;
Tmin = 0.001;
Tmax = 1;
low = Tmin;
upp = Tmax;
C_old = 0;


iter = 0;
while run
    iter = iter + 1;
    
    [K,F] = stiffOpt(U,V,INC,IEN,LM,gp,gw,p,q,B(:,1),B(:,2),T,B(:,3),D,R,dR_dx,J);
    F(ID(2,n*(m-1)+1)) = F(ID(2,n*(m-1)+1)) - 10000;
    
    d=solve(K,F);
    [dg0,dg1,vol] = Zsense(d,U,V,INC,IEN,LM,gp,gw,p,q,B(:,1),B(:,2),T,B(:,3),D,R,dR_dx,J);
    fval = vol-Vmax;
    C1 = d'*K'*d;
        [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp] = ... 
        mmasub(1,length(INC),iter,T,Tmin,Tmax,T_old1,T_old2, ...
        C1,dg0,fval,dg1,low,upp,1,0,1000,1);
    T_old2 = T_old1; 
    T_old1 = T;
    T = xmma;
    if norm(C1-C_old)<0.001 
        run = false;
    end
    C_old = C1;
end

%%

figure(24)
clf('reset')
hold on
nelm = (n-p)*(m-q);
X = zeros(steps,steps);
Y = zeros(steps,steps);
Z = zeros(steps,steps);
for e = 1:nelm
    A = IEN(e,:);
    Bex = B(A,1);
    Bey = B(A,2);
    Bez = T(A);
for j = 1:steps
    for k = 1:steps 
        [R, ~, ~] = Shape_function(pos(j), pos(k), e, p, q, B(:,1),B(:,2),B(:,3), U, V, INC, IEN);
        X(j,k) = Bex'*R;
        Y(j,k) = Bey'*R;
        Z(j,k) = Bez'*R;
        
    end
end
surf(X,Y,Z,'EdgeColor','none')
end
xlim([0,33]);
ylim([-2,12]);