kn = 101;
xu = linspace(0,4,kn);
p = 2; % Base function degree
U = [zeros(1,p+1),1,1,2,2,3,3,ones(1,p+1)*4]; % Open knot vector
n = length(U)-(p+1); % Number of base functions
Nu = zeros(n,kn); % Base function values
Bx = [2;2;0;-2;-2;-2;0;2;2]; %Control points
By = [0;2;2;2;0;-2;-2;-2;0];
b = 1/sqrt(2);
Bz = [1;b;1;b;1;b;1;b;1];
Bx = Bx.*Bz;
By = By.*Bz;

for j =1:kn %Calculates base functions values on the entire knot vector
    u = xu(j);
    i = FindSpan(n,p,u,U); 
    Nu(i+1-p:i+1,j) = BasisFuns(u,p,n,U);
end


figure(3)
clf('reset')
hold on

plot3(2*Bx,2*By,2*Bz,'ro-', 'MarkerFaceColor', 'g', 'MarkerEdgeColor' ,'g', 'Color' , 'k', 'LineWidth' , 1.5)
plot3(Bx./Bz,By./Bz,Bz./Bz, 'ro-', 'MarkerFaceColor', 'r', 'MarkerEdgeColor' ,'r', 'Color' , 'k', 'LineWidth' , 1.5)

X = zeros(3,kn);

for k = 1:kn %Coordinates for B-spline curve
    X(1,k) = Nu(:,k)'*2*Bx;
    X(2,k) = Nu(:,k)'*2*By;
    X(3,k) = Nu(:,k)'*2*Bz;
end

plot3(X(1,:),X(2,:),X(3,:), 'LineWidth' , 2);

for k = 1:10:kn
    xs = [0,X(1,k)];
    ys = [0,X(2,k)];
    zs = [0,X(3,k)];
    plot3(xs,ys,zs,'k')
end

xlim([-5,5]);
ylim([-5,5]);
r = 2;
x = 0;
y = 0;
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
zunit = ones(1,101);
plot3(xunit, yunit,zunit,'LineWidth',3.0);
surf([-5,5;-5,5],[-5,-5;5,5],[1,1;1,1]);
alpha 0.2