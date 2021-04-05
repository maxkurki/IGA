um = 5;
kn = um*20 + 1;
x = linspace(0,um,kn);
p = 2; % Base function degree
U = [zeros(1,p+1),1,2,3,3,4,ones(1,p+1)*um]; % Open knot vector
b = length(U)-(p+1); %Number of base functions
NN = zeros(b,kn); % Base function values
n = length(U) - p - 1; % Knot span index
for j =1:kn %Calculates base functions values on the entire knot vector
    u = x(j);
    i = FindSpan(n,p,u,U);
    N = BasisFuns2(u,p,n,0,U);
    NN(i+1-p:i+1,j) = N;
end
figure(1)
clf('reset');
hold on
for k = 1:b %Plots base functions
   plot(x,NN(k,:), 'LineWidth', 1.0);
end



Bx = [1;2;3;5;4;3;1;1];
By = [1;3;2;4;7;5;6;3];

[U, Bx, By, n, wn] = knotIns(U,2.5,p,n,Bx,By,1,ones(size(Bx)));

b = length(U)-(p+1); %Number of base functions
NN = zeros(b,kn); % Base function values
for j =1:kn %Calculates base functions values on the entire knot vector
    u = x(j);
    i = FindSpan(n,p,u,U);
    N = BasisFuns2(u,p,n,0,U);
    NN(i+1-p:i+1,j) = N;
end

figure(2)
clf('reset');
hold on
plot(Bx,By,'-bo','MarkerFaceColor','r','MarkerEdgeColor','r','LineWidth', 1.0);
xlim([0,6]);
ylim([0,8]);

X = zeros(2,kn);

for k = 1:kn %Coordinates for B-spline curve
    X(1,k) = NN(:,k)'*Bx;
    X(2,k) = NN(:,k)'*By;
end

plot(X(1,:),X(2,:),'-k','LineWidth',1.0) %Plot curve

figure(3)
clf('reset');
hold on
plot(X(1,:),X(2,:),'-k','LineWidth',1.0) %Plot curve
xlim([0,6]);
ylim([0,8]);
el = [1,1; 2.5,2.5; 4,3; 4.5,4.5; 4,7; 2,5.5; 1,3];
plot(el(:,1),el(:,2), 'rs', 'MarkerFaceColor', 'r');
