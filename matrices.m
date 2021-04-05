kn = 31;
x = linspace(0,1,kn);
p = 2; % Base function degree
q = 2;
U = [zeros(1,p+1),0.5,ones(1,p+1)*1]; % Open knot vector
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

%%

X = [0,0.5,1;0,0.5,1];
Y = [0,0,0;1,1,1];
Z = [1,1,1;1,1,1;];
figure(5)
clf('reset');
hold on
surf(X(:,1:2),Y(:,1:2),Z(:,1:2),'EdgeColor','none','FaceColor',[0.5 0 0.5]);
surf(X(:,2:3),Y(:,2:3),Z(:,2:3),'EdgeColor','none','FaceColor',[1 0 1]);
xticks([0 0.5 1]);
yticks([0 1]);
ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;

%%
figure(4)
clf('reset')
hold on

X = [5,6,6,5;3,4,4,3;4,5,5,4;5,6,6,5;3,4,4,3;4,5,5,4;5,6,6,5];
Y = [3,3,4,4;4,4,5,5;4,4,5,5;4,4,5,5;5,5,6,6;5,5,6,6;5,5,6,6];
fill(X',Y','r')
X = [1,2,2,1;2,3,3,2;3,4,4,3;1,2,2,1;2,3,3,2;3,4,4,3;1,2,2,1;2,3,3,2];
Y = [1,1,2,2;1,1,2,2;1,1,2,2;2,2,3,3;2,2,3,3;2,2,3,3;3,3,4,4;3,3,4,4];
fill(X',Y','b')
X = [3,4,4,3];
Y = [3,3,4,4];
fill(X',Y',[0.5, 0, 0.5])
X = [4,5,5,4];
Y = [3,3,4,4];
fill(X',Y',[1, 0, 1])
xlim([1,7]);
ylim([1,6]);
grid on
xticks([1 2 3 4 5 6 7])
xticklabels({'\xi_1','\xi_2','\xi_3','\xi_4','\xi_5','\xi_6','\xi_7'})
yticks([1 2 3 4 5 6])
yticklabels({'\eta_1','\eta_2','\eta_3','\eta_4','\eta_5','\eta_6'})
ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
%%
X = [-1;1;1;-1];
Y = [-1;-1;1;1];
figure(6)
clf('reset');
hold on
fill(X,Y,'b','Facealpha',0.6)
plot([-1.5,1.5],[0,0],'k')
plot([0,0],[-1.5,1.5],'k')
xlim([-1.5,1.5]);
xticks([-1 0 1])
yticks([-1 0 1])
ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
ylim([-1.5,1.5]);
grid on