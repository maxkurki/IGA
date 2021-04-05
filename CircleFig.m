kn = 73;
xu = linspace(0,4,kn);
p = 2; % Base function degree
U = [zeros(1,p+1),1,1,2,2,3,3,ones(1,p+1)*4]; % Open knot vector
xv = linspace(0,1,kn);
q = 1;
V = [zeros(1,q+1),ones(1,q+1)*1]; 
n = length(U)-(p+1); % Number of base functions
m = length(V)-(q+1); 
Nu = zeros(n,kn); % Base function values
Nv = zeros(m,kn);
for j =1:kn %Calculates base functions values on the entire knot vector
    u = xu(j);
    i = FindSpan(n,p,u,U); 
    Nu(i+1-p:i+1,j) = BasisFuns(u,p,n,U);
    v = xv(j);
    i = FindSpan(m,q,v,V); 
    Nv(i+1-q:i+1,j) = BasisFuns(v,q,m,V);
end

figure(3)
Bx = [2,1;2,1;0,0;-2,-1;-2,-1;-2,-1;0,0;2,1;2,1]; %Control points
By = [0,0;2,1;2,1;2,1;0,0;-2,-1;-2,-1;-2,-1;0,0];
b = 1/sqrt(2);
w = [1,1;b,b;1,1;b,b;1,1;b,b;1,1;b,b;1,1];


hold on
X = zeros(kn,kn);
Y = zeros(kn,kn);
Z = zeros(kn,kn);
W = NurbsW(Nu,Nv,w);

for j = 1:kn
for k = 1:kn %Coordinates for B-spline curve
    X(j,k) = Nu(:,j)'*(Bx.*w)*Nv(:,k)/W(j,k);
    Y(j,k) = Nu(:,j)'*(By.*w)*Nv(:,k)/W(j,k);
    Z(j,k) = sum(kron(Nu(:,j),Nv(:,j)));
end
end
figure(3)
hold on
surf(X,Y,Z,'EdgeColor','none','FaceColor','b');
xlim([-2.5,2.5]);
ylim([-2.5,2.5]);
r = 2;
x = 0;
y = 0;
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
plot(xunit, yunit,'k','LineWidth',3.0);
r = 1;
x = 0;
y = 0;
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
plot(xunit, yunit,'k','LineWidth',3.0);
for r = 1:length(Bx(:,1))
    plot3(Bx(r,:),By(r,:),ones(1,length(Bx(r,:))),'k-o','MarkerFaceColor','r','MarkerEdgeColor' , 'r', 'LineWidth', 1.0)
end
for r = 1:length(Bx(1,:))
    plot3(Bx(:,r),By(:,r),ones(1,length(Bx(:,r))),'k-o','MarkerFaceColor','r','MarkerEdgeColor' , 'r', 'LineWidth', 1.0)
end
plot(Bx,By,'ro','MarkerFaceColor','r');