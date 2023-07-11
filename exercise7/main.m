A = [5, 1, 1;
    1, 5, -1;
    1, -1, 5];
b = [4 2 -4]';
x0 = [0 0 0]';

%% a) ii
A1 = matrix_7_3(100,10);
A2 = matrix_7_3(1000,100);
b1 = ones(100,1);
b2 = ones(1000,1);
x01 = zeros(100,1);
x02 = zeros(1000,1);
tol = 1e-10;
maxit = 5000;

[x1,k1,resvec1] = cg(A1,b1,x01,tol,maxit);
[x2,k2,resvec2] = cg(A2,b2,x02,tol,maxit);
disp("Matrix of size (n,m) = (100,10) converges after k=" + k1 + " iterations.")
disp("Matrix of size (n,m) = (1000,100) converges after k=" + k2 + " iterations.")

%% c)
loglog(1:k1, resvec1)
title('CG convergence study for $(n,m) = (100,10)$', 'Interpreter', 'Latex')
xlabel('Iteration number $k$', 'Interpreter', 'Latex')
ylabel('Absolute $L^2$ error $||\mathbf{Ax}^{(k)} - \mathbf{b}||_2$', 'Interpreter', 'Latex')
legend("Absolute residual error")
grid on

figure
loglog(1:k2, resvec2)
title('CG convergence study for $(n,m) = (1000,100)$', 'Interpreter', 'Latex')
xlabel('Iteration number $k$', 'Interpreter', 'Latex')
ylabel('Absolute $L^2$ error $||\mathbf{Ax}^{(k)} - \mathbf{b}||_2$', 'Interpreter', 'Latex')
legend("Absolute residual error")
grid on