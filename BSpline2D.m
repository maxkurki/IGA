%Example of B-spline surface

kn = 31;
x = linspace(0,1,kn);
p = 2; % Base function degree
q = 2;
U = [zeros(1,p+1),0.5,0.5,ones(1,p+1)*1]; % Open knot vector
V = [zeros(1,q+1),ones(1,q+1)*1]; 
bu = length(U)-(p+1); % Number of base functions
bv = length(V)-(q+1); 
Nu = zeros(bu,kn); % Base function values
Nv = zeros(bv,kn);
n = length(U) - p - 1; % Knot span index
m = length(V) - q - 1;
for j =1:kn %Calculates base functions values on the entire knot vector
    u = x(j);
    i = FindSpan(n,p,u,U); 
    Nu(i+1-p:i+1,j) = BasisFuns2(u,p,n,0,U);
    i = FindSpan(m,q,u,V); 
    Nv(i+1-p:i+1,j) = BasisFuns2(u,q,m,0,V);
end

figure(1)
clf('reset');
hold on
for k = 1:bu %Plots base functions
   plot(x,Nu(k,:));
end

figure(2)
clf('reset');
hold on
for k = 1:bv %Plots base functions
   plot(-x,Nv(k,:));
end

figure(3)
clf('reset');
hold on
Bx = [0,-1,-2;0,-1,-2;0.5,-1,-2;1,1,1;2,2,2]; %Control points
By = [0,0,0;1,1,1;1.5,3,4;2,3,4;2,3,4];
for r = 1:5
    plot3(Bx(r,:),By(r,:),ones(1,length(Bx(r,:))),'k-o','MarkerFaceColor','r','MarkerEdgeColor' , 'r', 'LineWidth', 1.0)
end
for r = 1:3
    plot3(Bx(:,r),By(:,r),ones(1,length(Bx(:,r))),'k-o','MarkerFaceColor','r','MarkerEdgeColor' , 'r', 'LineWidth', 1.0)
end
%plot(Bx,By,'ro','MarkerFaceColor','r');
xlim([-2.5,2.5]);
ylim([-1,5]);

X = zeros(kn,kn);
Y = zeros(kn,kn);
Z = zeros(kn,kn);

for j = 1:kn
for k = 1:kn %Coordinates for B-spline curve
    X(j,k) = Nu(:,j)'*Bx*Nv(:,k);
    Y(j,k) = Nu(:,j)'*By*Nv(:,k);
    Z(j,k) = sum(kron(Nu(:,j),Nv(:,j)));
end
end

figure(3)
%clf('reset');
surf(X(1:ceil(kn/2),:),Y(1:ceil(kn/2),:),Z(1:ceil(kn/2),:),'EdgeColor','none','FaceColor','b')
hold on
surf(X(ceil(kn/2):end,:),Y(ceil(kn/2):end,:),Z(ceil(kn/2):end,:),'EdgeColor','none','FaceColor','r')
xlim([-2.5,2.5]);
ylim([-1,5]);
grid off

X = [0,0.5,1;0,0.5,1];
Y = [0,0,0;1,1,1];
Z = [1,1,1;1,1,1;];
figure(5)
clf('reset');
hold on
surf(X(:,1:2),Y(:,1:2),Z(:,1:2),'EdgeColor','none','FaceColor','b');
surf(X(:,2:3),Y(:,2:3),Z(:,2:3),'EdgeColor','none','FaceColor','r');
