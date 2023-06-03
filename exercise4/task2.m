%% Exercise b)
kmax = 100;
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
error_none = zeros(kmax, 1);
error_naive = zeros(kmax, 1);
error_wilkinson = zeros(kmax, 1);
for k = 1:kmax
    error_none(k) = norm(sort(ev) - sort(diag(iter_none(:,:,k))), 1);
    error_naive(k) = norm(sort(ev) - sort(diag(iter_naive(:,:,k))), 1);
    error_wilkinson(k) = norm(sort(ev) - sort(diag(iter_wilkinson(:,:,k))), 1);
end

semilogy(1:kmax, error_none, "g-")
hold on
grid on
semilogy(1:kmax, error_naive, "r--")
semilogy(1:kmax, error_wilkinson, "b-.")
legend("No shift", "Naive shift", "Wilkinson shift")
xlabel('Iteration number $k$', 'Interpreter', 'Latex')
ylabel('Absolute error $e_k = || \lambda - \lambda^k ||_1$', 'Interpreter', 'Latex')