% Devin Balian 2791430
clear all;

A = [6, 2, 1;
    2, 3, 1;
    1, 1, 1];
mu = 3;
x0 = [1;1;1];
kmax = 15;

[lambda, x] = inv_veciter(A, mu, x0, kmax)

[V, D] = eig(A)