% Devin Balian 2791430
clear all;

%%% Exercise a) iii)
A = [1, 2;
    2, 3];
kmax = 2;

M = magic(5)
G = givens(M, 1, 2, 2)
G * M

ev_num = qr_algorithm(A, 'wilkinson', kmax);
ev_exact = eig(A);

ev_error = norm(sort(ev_num) - sort(ev_exact)) / norm(ev_exact);
disp(ev_error)