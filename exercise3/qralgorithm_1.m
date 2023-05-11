% Devin Balian 2791430

function [ev] = qralgorithm_1(A, kmax)

    % QR Iteration
    for k = 1:kmax
        [Q, R]= qr(A);
        A = R*Q;
    end

    % Get diagonal elements of A
    ev = diag(A);
end