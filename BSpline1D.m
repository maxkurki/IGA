%Example of B-spline curve

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



B = [1,1;2,3;3,2;5,4;4,7;3,5;1,6;1,3];
figure(2)
clf('reset');
hold on
plot(B(:,1),B(:,2),'-bo','MarkerFaceColor','r','MarkerEdgeColor','r','LineWidth', 1.0);
xlim([0,6]);
ylim([0,8]);


for k = 1:kn %Coordinates for B-spline curve
    X(1,k) = NN(:,k)'*B(:,1);
    X(2,k) = NN(:,k)'*B(:,2);
end

plot(X(1,:),X(2,:),'-k','LineWidth',1.0) %Plot curve

figure(3)
clf('reset');
hold on
plot(X(1,:),X(2,:),'-k','LineWidth',1.0) %Plot curve
xlim([0,6]);
ylim([0,8]);
plot(X(1,1:floor(kn/um):end),X(2,1:floor(kn/um):end), 'rs', 'MarkerFaceColor', 'r');
