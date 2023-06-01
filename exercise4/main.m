% Devin Balian 2791430
clear all;

%% Exercise a) iii)
A = [1, 2;
    2, 3];
kmax = 2;

ev_none = qr_algorithm(A, 'none', kmax);
ev_naive = qr_algorithm(A, 'naive', kmax);
ev_wilkinson = qr_algorithm(A, 'wilkinson', kmax);
ev = eig(A);

error_none = norm(sort(ev_none) - sort(ev)) / norm(ev);
error_naive = norm(sort(ev_naive) - sort(ev)) / norm(ev);
error_wilkinson = norm(sort(ev_wilkinson) - sort(ev)) / norm(ev);

disp("Error without shift: " + error_none)
disp("Error with naive shift: " + error_naive)
disp("Error with wilkinson shift: " + error_wilkinson)

%% Exercise b)
kmax = 200;
N = 50;
D = diag(1:N);
[Q,~] = qr(2*rand(N,N) - ones(N,N));
A = Q*D*Q';

ev = 1:N;
tic
[~, iter_none] = qr_algorithm(A, 'none', kmax);
toc
[~, iter_naive] = qr_algorithm(A, 'naive', kmax);
[~, iter_wilkinson] = qr_algorithm(A, 'wilkinson', kmax);

% Calculate errors
error_none = zeros(kmax, 1);
error_naive = zeros(kmax, 1);
error_wilkinson = zeros(kmax, 1);
for k = 1:kmax
    error_none(k) = norm(sort(ev) - sort(diag(iter_none(:,:,k)))', 1);
    error_naive(k) = norm(sort(ev) - sort(diag(iter_naive(:,:,k)))', 1);
    error_wilkinson(k) = norm(sort(ev) - sort(diag(iter_wilkinson(:,:,k)))', 1);
end

semilogy(1:kmax, error_none, "g-")
hold on
grid on
semilogy(1:kmax, error_naive, "r--")
semilogy(1:kmax, error_wilkinson, "b-.")
legend("No shift", "Naive shift", "Wilkinson shift")
xlabel('Iteration number $k$', 'Interpreter', 'Latex')
ylabel('Absolute error $e_k = || \lambda - \lambda^k ||_1$', 'Interpreter', 'Latex')