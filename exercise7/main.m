% Devin Balian 2791430
clear all

%% a) ii
A1 = matrix_7_3(100,10);
A2 = matrix_7_3(1000,100);
b1 = ones(100,1);
b2 = ones(1000,1);
x01 = zeros(100,1);
x02 = zeros(1000,1);
tol = 1e-10;
maxit = 5000;

[~,k1,resvec1] = cg(A1,b1,x01,tol,maxit);
[~,k2,resvec2] = cg(A2,b2,x02,tol,maxit);
disp("task a)")
disp("Matrix of size (n,m) = (100,10) converges after k=" + k1 + " iterations.")
disp("Matrix of size (n,m) = (1000,100) converges after k=" + k2 + " iterations.")

%% b) ii
P1 = diag(diag(A1));
P2 = diag(diag(A2));
[~,k1p,resvec1p] = pcg(A1,b1,P1,x01,tol,maxit);
[~,k2p,resvec2p] = pcg(A2,b2,P2,x02,tol,maxit);
disp("task b)")
disp("Preconditioned matrix of size (n,m) = (100,10) converges after k=" + k1p + " iterations.")
disp("Preconditioned matrix of size (n,m) = (1000,100) converges after k=" + k2p + " iterations.")

%% c)
c1 = condest(A1);
c1p = condest(inv(P1)*A1);
c2 = condest(A2);
c2p = condest(inv(P2)*A2);

disp("task c)")
disp("Condition of matrix of size (n,m) = (100,10):")
disp("  cond(A) = " + c1);
disp("  cond(P^-1 A) = " + c1p);
disp("Condition of matrix of size (n,m) = (1000,100):")
disp("  cond(A) = " + c2);
disp("  cond(P^-1 A) = " + c2p);

loglog(1:k1, resvec1)
hold on
loglog(1:k1p, resvec1p)
title('CG convergence study for $(n,m) = (100,10)$', 'Interpreter', 'Latex')
xlabel('Iteration number $k$', 'Interpreter', 'Latex')
ylabel('Absolute $L^2$ error $||\mathbf{Ax}^{(k)} - \mathbf{b}||_2$', 'Interpreter', 'Latex')
yline(tol, "r--")
legend("CG", "PCG", "Error tolerance", 'Location', 'west')
grid on

figure
loglog(1:k2, resvec2)
hold on
loglog(1:k2p, resvec2p)
title('CG convergence study for $(n,m) = (1000,100)$', 'Interpreter', 'Latex')
xlabel('Iteration number $k$', 'Interpreter', 'Latex')
ylabel('Absolute $L^2$ error $||\mathbf{Ax}^{(k)} - \mathbf{b}||_2$', 'Interpreter', 'Latex')
yline(tol, "r--")
legend("CG", "PCG", "Error tolerance", 'Location', 'west')
grid on