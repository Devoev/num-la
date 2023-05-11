% Devin Balian 2791430
clear all;

%%% Exercise a) ii)
A = [10, -6, 8;
    -6, 15, 10;
    8, 10, 5];
kmax = 2000;
tol = 1e-14;

% Calculate eigenvalues and relative errors
ev1 = qralgorithm_1(A, kmax);
[ev2, iter] = qralgorithm_2(A, tol, kmax);
ev = eig(A);

e1 = norm(sort(ev1) - sort(ev)) / norm(ev);
e2 = norm(sort(ev2) - sort(ev)) / norm(ev);

disp("a) ii)")
disp("Relative error of ev1: " + e1 + ". Stopped after " + kmax + " iterations.")
disp("Relative error of ev2: " + e2 + ". Converged after " + iter + " iterations.")

%%% Exercise a) iii)
A = @(a) [a, 2, 3, 13;
        5, 11, 10, 8;
        9, 7, 6, 12;
        4, 14, 15, 1];

ev1 = eig(A(30));
[evQR1, iter1] = qralgorithm_2(A(30), tol, kmax);
ev2 = eig(A(-30));
[evQR2, iter2] = qralgorithm_2(A(-30), tol, kmax);

e1 = norm(sort(evQR1) - sort(ev1)) / norm(ev1);
e2 = norm(sort(evQR2) - sort(ev2)) / norm(ev2);

disp("a) iii)")
disp("Relative error of ev1 (A(30)): " + e1 + ". Converged after " + iter1 + " iterations.")
disp("Relative error of ev2 (A(-30)): " + e2 + ". Converged after " + iter2 + " iterations.")