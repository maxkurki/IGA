%Plots B-spline basis functions of order 0,1,2

figure(1)
clf('reset');
hold on
plot([0,1],[1,1],'k--', 'LineWidth', 1.5)
plot([1,2],[1,1],'r--', 'LineWidth', 1.5)
plot([2,3],[1,1],'b--', 'LineWidth', 1.5)
plot([0,5],[0,0],'k-')
xlim([0,5])
ylim([0,4])
xticks([0 1 2 3 4 5])
yticks([0 1])
%%
figure(2)
clf('reset');
hold on
plot([0 1 2],[0 1 0],'k--', 'LineWidth', 1.5);
plot([1 2 3],[0 1 0],'r--', 'LineWidth', 1.5);
plot([2 3 4],[0 1 0],'b--', 'LineWidth', 1.5);
plot([0,5],[0,0],'k-')
xlim([0,5])
ylim([0,4])
xticks([0 1 2 3 4 5])
yticks([0 1])
%%
figure(3)
clf('reset');
hold on
x = linspace(0,0.99,100);
y = (x.^2)/2;
x = linspace(1,1.99,100);
y = [y, (-2*x.^2 + 6*x - 3)/2];
x = linspace(2,2.99,100);
y = [y, ((3 - x).^2)/2];
plot(linspace(0,2.99,300),y,'k--', 'LineWidth', 1.5);
plot(linspace(1,3.99,300),y,'r--', 'LineWidth', 1.5);
plot(linspace(2,4.99,300),y,'b--', 'LineWidth', 1.5);
plot([0,5],[0,0],'k-')
xlim([0,5])
ylim([0,4])
xticks([0 1 2 3 4 5])
yticks([0 1])
