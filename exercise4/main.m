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