function [V,H] = arnoldi(A,r0,m,tol)
% ARNOLDI Calculates an orthonomral basis of the Krylov space K^m(A,r0).
% Inputs:
%   A   - System matrix of size (n,n).
%   r0  - Initial basis vector of size (n).
%   m   - Dimension of the Krylov space.
%   tol - Tolerance.
% Outputs:
%   V   - Orthonormal basis of size (n,m+1).
%   H   - Hessenberg matrix of size (m+1,m).

    % Init matrices
    [n,~] = size(A);
    V = zeros(n,m+1);
    V(:,1) = r0 / norm(r0);
    H = zeros(m+1,m);

    for j = 1:m
        % Arnoldi step
        z = A*V(:,j);

        % Orthonormalization
        for i = 1:j
            H(i,j) = V(:,i)' * z;
            z = z - H(i,j)*V(:,i);
        end

        H(j+1,j) = norm(z);

        % Check for stagnation
        if abs(H(j+1,j)) < tol
            return
        else
            V(:,j+1) = z / H(j+1,j);
        end
    end

end

% Devin Balian 2791430