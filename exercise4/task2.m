%% Exercise b)
kmax = 500;
tol = 1e-10;
N = 100;
D = diag(1:N);
[Q,~] = qr(2*rand(N,N) - ones(N,N));
A = Q*D*Q';
%A = full(gallery('tridiag',N,-1,2,-1));
%A = wilkinson(N, 'single');

ev = eig(A);
[~, iter_none] = qr_algorithm(A, 'none', kmax, tol, false);
[~, iter_naive] = qr_algorithm(A, 'naive', kmax, tol, false);
[~, iter_wilkinson] = qr_algorithm(A, 'wilkinson', kmax, tol, false);

% Calculate errors
error_none = vecnorm(sort(ev) - sort(iter_none), 1);
error_naive = vecnorm(sort(ev) - sort(iter_naive), 1);
error_wilkinson = vecnorm(sort(ev) - sort(iter_wilkinson), 1);

semilogy(1:kmax, error_none, "g-")
hold on
grid on
semilogy(1:kmax, error_naive, "r--")
semilogy(1:kmax, error_wilkinson, "b-.")
legend("No shift", "Naive shift", "Wilkinson shift")
xlabel('Iteration number $k$', 'Interpreter', 'Latex')
ylabel('Absolute error $e_k = || \lambda - \lambda^k ||_1$', 'Interpreter', 'Latex')