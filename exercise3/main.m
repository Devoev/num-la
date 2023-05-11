% Devin Balian 2791430
clear all;

A = [10, -6, 8;
    -6, 15, 10;
    8, 10, 5];
kmax = 2000;
tol = 1e-14;

ev1 = qralgorithm_1(A, kmax)
ev2 = qralgorithm_2(A, tol, kmax)
ev = eig(A)