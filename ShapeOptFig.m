%Shape optimization of a cantilever

E = 2.1*10^11; % Young's modulus in N/m2
v = 0.3; % Poission's ratio

gp = [-sqrt(3/5),0,sqrt(3/5)]; %Gauss quadrature points
gw = [5/9,8/9,5/9]; %Guass quadrature point weights

D = hooke(1,E,v); % Constitutive matrix for plane stress

p = 2; % Base function degrees
q = 2;
U = [zeros(1,p+1),1,2,3,ones(1,p+1)*4]; % Open knot vectors
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

Bx = kron(linspace(0,30,n)',ones(1,m)); %Uniform control point grid as nxm matrix
By = kron(ones(n,1),linspace(0,10,m));

Bc = zeros(n*m,3);
for i = 1:m %Control points as vector and control point weights
    for j = 1:n
        Bc(j+(i-1)*n,1) = Bx(j,i);
        Bc(j+(i-1)*n,2) = By(j,i);
        Bc(j+(i-1)*n,3) = 1;
    end
end

[C,Cp] = CP(INC,IEN,m); %Connectivity matrix for shape optimization

mx = 8; %Allowed movement in x- and y-direction respectively
my = 8.5;
Bc(Cp(:,1),2) = Bc(Cp(:,1),2) - my; % "Constants" of control points

Vmax = 0.7*10*30; % Maximum allowed volume


T = 0.5*ones(max(C(:,1)),1); %Optimization parameters

B = Bc; %Control points values
B(Cp(:,1),2) = B(Cp(:,1),2) + T*my;
B(Cp(:,1)-n,2) = B(Cp(:,1),2)/2;

B(15,2) = B(15,2) + 2;

T_old1 = T; %MMA parameters
T_old2 = T;
Tmin = 0.001;
Tmax = 1;
low = Tmin;
upp = Tmax;
C_old = 0;

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


iter = 0;
run = true;
while run
    iter = iter + 1;
    
    [K,F] = stiff(U,V,INC,IEN,LM,gp,gw,p,q,B(:,1),B(:,2),B(:,3),D); %Stiffness matrix and force vector
    F(ID(2,n)) = F(ID(2,n)) - 10000;
    d=solve(K,F); %Solve for displacements
    [dg0,dg1,vol] = Bsense(d,U,V,INC,IEN,LM,gp,gw,p,q,B(:,1),B(:,2),B(:,3),D,C,mx,my); %Sensitivities and volume
    fval = vol-Vmax; %current value of constraint function
    C1 = d'*K'*d; %Compliance
    [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp] = ... %Use MMA to clalculate new values
        mmasub(1,length(Cp(:,1)),iter,T,Tmin,Tmax,T_old1,T_old2, ...
        C1,dg0,fval,dg1,low,upp,1,0,1000,1);
    T_old2 = T_old1; %Update values
    T_old1 = T;
    T = xmma;
    B = Bc; %Update control points
    B(Cp(:,1),2) = B(Cp(:,1),2) + T*my;
    B(Cp(:,1)-n,2) = B(Cp(:,1),2)/2;
    if norm(C1-C_old)<0.0001 %Check convergence
        run = false;
    end
    C_old = C1;
end

%%

figure(25) %Plot optimized structure
clf('reset')
hold on
nelm = (n-p)*(m-q); %number of elements
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
        [R, ~, ~] = Shape_function(pos(j), pos(k), e, p, q, B(:,1),B(:,2),B(:,3), U, V, INC, IEN); %Basis function values in current evaluation point
        X(j,k) = Bex'*R; %Coordinates resulting from current evaluation point
        Y(j,k) = Bey'*R;
        
    end
end
surf(X,Y,Z,'EdgeColor','none') %Plot surface
end
xlim([0,33]);
ylim([-2,12]);

plot(B(:,1),B(:,2),'ro','MarkerFaceColor','r','MarkerEdgeColor', 'k'); %Plot control points
