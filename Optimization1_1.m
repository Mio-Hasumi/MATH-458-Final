% code for plotting the 2‑simplex {x1>=0, x2>=0, x1+x2<=1}

% Vertices of the triangle in (x1,x2):
V = [0 0; 1 0; 0 1];

figure; 
fill(V(:,1), V(:,2), 0.8*[1 1 1], 'LineWidth',1.5);
hold on

% Draw the constraint lines:
plot([0 1],[0 0],'k','LineWidth',1);     % x2 = 0
plot([0 0],[0 1],'k','LineWidth',1);     % x1 = 0
plot([1 0],[0 1],'k','LineWidth',1);     % x1 + x2 = 1

% in (x1,x2) we have z = 3 − 2x1 − x2.
c_vals = [1 1.5 2 2.5 3];  % pick some cost‐levels
x1 = linspace(0,1,200);
x2 = linspace(0,1,200);
[X1,X2] = meshgrid(x1,x2);
Z = 3 - 2*X1 - 1*X2;
contour(X1,X2, Z, c_vals, 'LineColor','b','LineStyle','--');

% Mark the three corner solutions:
pts = [1 0; 0 1; 0 0];
scatter(pts(:,1),pts(:,2),80,'ro','filled')
text(pts(:,1)+0.02, pts(:,2)+0.02, {'(1,0,0)','(0,1,0)','(0,0,1)'},...
     'FontSize',12)

axis equal tight

title('Feasible region: $x_1\ge0,\;x_2\ge0,\;x_1+x_2\le1$','Interpreter','latex')
xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
legend({'Feasible simplex','Constraints','Objective level sets','Vertices'}, ...
       'Interpreter','latex','Location','northeast')
grid on
