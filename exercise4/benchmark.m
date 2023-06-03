%% Benchmark matrix of exercise b)
nmax = 10; % Maximum number of times, the algorithm is executed
kmax = 200; % Maximum number of iterations of each qr algorithm
tol = 1e-10;
N = 100;
D = diag(1:N);
[Q,~] = qr(2*rand(N,N) - ones(N,N));
A = Q*D*Q';

tic
for n = 1:nmax
    qr_algorithm(A, 'none', kmax, tol, true);
end
toc