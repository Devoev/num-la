% Devin Balian 2791430
clear all;

%%% Exercise a) iii)
A = [1, 2;
    2, 3];
kmax = 2000;

ev = qralgorithm(A, 'abcd', kmax)