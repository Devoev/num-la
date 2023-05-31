% Devin Balian 2791430

function [ev, iter] = qralgorithm_1(A, tol, kmax)

    % QR Iteration
    iter = kmax;
    for k = 1:kmax
        % Get max absolute value of lower triangular part and check tolerance
        if max(abs(tril(A,-1)), [], 'all') <= tol
            iter = k;
            break
        end
        [Q, R]= qr(A);
        A = R*Q;
    end

    % Get diagonal elements of A
    ev = diag(A);
end