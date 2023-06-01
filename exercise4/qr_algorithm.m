function [ev] = qr_algorithm(A, shift, kmax)
% QR_ALGORITHM Caluclates the eigenvalues of A using the QR algorithm.
% Inputs:
%   shift - The shift technique. Use 'none', 'naive' or 'wilkinson'.
%   kmax  - The maximum amount of iterations.

    % Transformation to hessenberg form
    H = hess(A);
    [n,~] = size(H);

    % QR Iteration
    for k = 1:kmax
        % Calculate shift sigma
        if strcmp(shift, 'none')
            sigma = 0;
        elseif strcmp(shift, 'naive')
            sigma = H(n, n);
        elseif strcmp(shift, 'wilkinson')
            d = (H(n-1,n-1) - H(n,n)) / 2;
            sigma = H(n,n) + d - sign(d)*sqrt(d^2 + H(n-1,n)^2);
        else
            error("Unsupported shift technique. Use 'none', 'naive' or 'wilkinson'.");
        end

        [Q, R]= qr(H - sigma*eye(n));
        H = R*Q + sigma*eye(n);
    end

    % Get diagonal elements of A
    ev = diag(H);
end

% Devin Balian 2791430