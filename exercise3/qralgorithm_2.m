% Devin Balian 2791430

function [ev, iter] = qralgorithm_1(A, tol, kmax)

    % QR Iteration
    for k = 1:kmax
        if max(abs(tril(A,-1)), [], 'all') <= tol
            disp("Tolerance reached")
            break
        end
        [Q, R]= qr(A);
        A = R*Q;
    end

    % Get diagonal elements of A
    ev = diag(A);
end