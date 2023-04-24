% Devin Balian 2791430

clear all;

N = 15; % number of iterations
delta = 1e-5; % maximum disturbance error
e = zeros(N, 1); % ratio of e_x and e_b
cond = zeros(N, 1); % condition number of hilbert matrix

for n = 1:N
    % Hilbert matrix
    H = hilb(n);
    Hinv = invhilb(n);

    % rhs vector
    b = ones(n, 1);
    b_delta = b + rand(n, 1) * delta / sqrt(n);
    e_b = norm(b - b_delta)/norm(b);

    % Solve linear system
    x = H\b;
    x_delta = H\b_delta;
    e_x = norm(x - x_delta)/norm(x);

    % Set error and cond
    e(n) = e_x/e_b;
    cond(n) = norm(H)*norm(Hinv);
end

% Calculate mean value of a (=alpha) for cond(n) = exp(an)*exp(-1) = exp(an-a)
a = mean(diff(log(cond)));
c = exp((1:N) * a);
% c = c / exp(a); % optional y-axis shift

% Plot resuts
semilogy(e)
hold on
semilogy(cond)
semilogy(c)
legend("e_x/e_b", "cond(H)", "c(n)")
