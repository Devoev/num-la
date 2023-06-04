function [ev, iter] = qr_algorithm(A, shift, kmax, tol, deflation)
% QR_ALGORITHM Caluclates the eigenvalues of A using the QR algorithm.
% Inputs:
%   shift     - The shift technique. Use 'none', 'naive' or 'wilkinson'.
%   kmax      - The maximum amount of iterations.
%   tol       - Tolerance value.
%   deflation - Whether to apply deflation (only works on real symmetric matrices).
% Outputs:
%   ev        - Array of approximated eigenvalues of A.
%   iter      - The iteration of apprxoimated eigenvalues.

    % Test if A is a scalar
    if isscalar(A)
        ev = A;
        return
    end

    % Transformation to hessenberg form
    [n,~] = size(A);
    H = hess(A);
%    H = spdiags(spdiags(H, -1:1), -1:1, n, n);
%    H = full(H);
    iter = zeros(n, kmax);

    % QR Iteration
    for k = 1:kmax
        % Calculate shift s
        if strcmp(shift, 'none')
            s = 0;
        elseif strcmp(shift, 'naive')
            s = H(n,n);
        elseif strcmp(shift, 'wilkinson')
            d = (H(n-1,n-1) - H(n,n)) / 2;
            s = H(n,n) + d - sign(d)*sqrt(d^2 + H(n-1,n)^2);
        else
            error("Unsupported shift technique. Use 'none', 'naive' or 'wilkinson'.");
        end

        % Deflation
        if deflation
            lower_diag = diag(H, -1);
            idx = abs(lower_diag) < tol;
            idx = nonzeros(double(idx)' .* (1:n-1)); % Find all zero indices on lower diagonal
            if ~isempty(idx)
                i = idx(ceil(end/2)); % Get most centered zero index
                if i == 1 | i == n
                    break
                end
                disp("Calling deflation on index " + i)
                D1 = H(1:i-1, 1:i-1);
                D2 = H(i:n, i:n);
                ev = [qr_algorithm(D1, shift, kmax-k, tol, deflation); qr_algorithm(D2, shift, kmax-k, tol, deflation)];
                return
            end
        end

%        [Q, R]= qr_givens(H - s*eye(n), tol, deflation);
        [Q, R]= qr(H - s*eye(n));
        H = R*Q + s*eye(n);
        iter(:,k) = diag(H);
    end

    % Get eigenvalues elements of last iteration
    ev = diag(H);
end

% Devin Balian 2791430