% Devin Balian 2791430

clear all;

N = 15;
delta = 1e-5;
e = zeros(N, 1);
cond = zeros(N, 1);
for n = 1:N

    % Hilbert matrix
    H = hilb(n);
    Hinv = invhilb(n);

    % rhs vector
    b = ones(n, 1);
    b_delta = b + rand(n, 1) * delta / sqrt(n);
    e_b = norm(b - b_delta)/norm(b);

    x = H\b;
    x_delta = H\b_delta;
    e_x = norm(x - x_delta)/norm(x);

    % Set error and data
    e(n) = e_x/e_b;
    cond(n) = norm(H)*norm(Hinv);
end

% Calculate a for cond(n) = exp(an)
a = mean(diff(log(cond)))
c = exp((1:N) * a) / exp(a);

% Plot resuts
semilogy(e)
hold on
semilogy(cond)
semilogy(c)
legend("e_x/e_b", "cond(H)", "c(n)")
