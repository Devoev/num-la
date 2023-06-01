% Devin Balian 2791430

function [ev] = qralgorithm(A, shift, kmax)
    % Transformation to hessenberg form
    H = hess(A);

    % QR Iteration
    for k = 1:kmax
        % Calculate shift sigma
        if strcmp(shift, 'none')
            sigma = 0;
        elseif strcmp(shift, 'naive')
            sigma = H(end, end);
        elseif strcmp(shift, 'wilkinson')
            sigma = None;
        else
            error("Unsupported shift routine. Use 'none', 'naive' or 'wilkinson'.");
        end

        [Q, R]= qr(H);
        H = R*Q;
    end

    % Get diagonal elements of A
    ev = diag(H);
end