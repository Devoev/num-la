% Devin Balian 2791430
clear all;

A = [6, 2, 1;
    2, 3, 1;
    1, 1, 1];
mu = 3;
x0 = [1;1;1];
kmax = 15;

% Inverse vector iteration
[lambda, x] = inv_veciter(A, mu, x0, kmax);

% "exact" eigenvalues/vectors
[V, D] = eig(A);
lambda3 = D(1,1);
lambda2 = D(2,2);
lambda1 = D(3,3);
x2 = V(:,2);

% Calculate errors
k_range = 1:kmax;
e_lambda = abs(lambda2 - lambda);
e_vec = arrayfun(@(k) min( norm(x2 - x(:,k), Inf) , norm(x2 + x(:,k), Inf) ), k_range);

% Additional lines
c1 = power(abs((lambda2 - mu) / (lambda3 - mu)), 2*k_range);
c2 = power(abs((lambda2 - mu) / (lambda1 - mu)), 2*k_range);
c3 = power(abs((lambda2 - mu) / (lambda3 - mu)), k_range);
c4 = power(abs((lambda2 - mu) / (lambda1 - mu)), k_range);

% Plot errors
semilogy(k_range, e_lambda,...
    k_range, e_vec,...
    k_range, c1, "--",...
    k_range, c2, "--",...
    k_range, c3, "--",...
    k_range, c4, "--")
xlabel("Iteration k")
title("Inverse vector iteration convergence analysis")
legend("Absolute error (2nd eigenvalue)", "Inf norm error (2nd eigenvector)", "c1(k)", "c2(k)", "c3(k)", "c4(k)")
grid on