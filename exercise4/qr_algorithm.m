function [ev, iter] = qr_algorithm(A, shift, kmax, tol, varargin)
% QR_ALGORITHM Caluclates the eigenvalues of A using the QR algorithm.
% Inputs:
%   shift     - The shift technique. Use 'none', 'naive' or 'wilkinson'.
%   kmax      - The maximum amount of iterations.
%   tol       - Tolerance value.
% Varargs
%   deflation - Apply deflation (only works on real symmetric matrices).
%   tridiag   - Transform A to hessenberg form (tridiagonal in the symmetric case).
% Outputs:
%   ev        - Array of approximated eigenvalues of A.
%   iter      - The iteration of apprxoimated eigenvalues.

    % Test if A is a scalar
    if isscalar(A)
        ev = A;
        iter = A;
        return
    end

    % Get options
    option_deflate = any(strcmp(varargin,'deflation'));
    option_tridiag = any(strcmp(varargin,'tridiag'));

    % Transformation to hessenberg form
    if option_tridiag
        H = hess(A);
    else
        H = A;
    end

    [n,~] = size(A);
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

        % Deflation - test if tridiagonal element is close to zero
        if option_deflate & abs(H(n, n-1)) < tol * min(abs(H(n-1, n-1)), abs(H(n,n)))
            D1 = H(1:n-1, 1:n-1);
            ev2 = H(n,n);
            [ev1, iter1] = qr_algorithm(D1, shift, kmax-k+1, tol, varargin{:}); % Calculate eigenvalues of upper block

            % Merge iterated eigenvalues
            iter(1:end-1, k:end) = iter1;
            iter(end, k:end) = ev2;
            ev = [ev1; ev2];
            return
        end

%        [Q, R]= qr_givens(H - s*eye(n), tol, option_tridiag);
        [Q, R]= qr(H - s*eye(n)); % Using internal Matalb QR decomposition, because qr_givens takes way too long for large matrices
        H = R*Q + s*eye(n);
        iter(:,k) = diag(H);
    end

    % Get eigenvalues elements of last iteration
    ev = diag(H);
end

% Devin Balian 2791430