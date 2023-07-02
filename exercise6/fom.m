function [x,iter] = fom(A,b,kmax,tol)
% FOM Calculates an approximate solution of Ax=b using the Full Orthogonalization Method (FOM).
% Inputs:
%   A    - System matrix of size (n,n).
%   b    - Right hand side of size(n).
%   kmax - Maximum number of iterations.
%   tol  - Tolerance.
% Outputs:
%   x    - Approximate solution of size (n).
%   iter - Iteration of solutions.

    [n,~] = size(A);
    iter = zeros(n,kmax);
    for k=2:kmax
        [V,H] = arnoldi(A,b,k,tol);
        Vk = V(:,1:k);
        Hk = H(1:k,:);

        % Solve H-system
        rhs = norm(b) * eye(k,1);
        y = linsolve(Hk,rhs);

        % Backwards Projection
        iter(:,k) = Vk*y;
        x = iter(:,k);

        % Termination criteria
        if norm(A*x - b) < tol | norm(x - iter(:,k-1)) < tol
            iter(:,k:kmax) = repmat(x,1,kmax-k+1);
            disp("FOM converged after k=" + k + " iterations.")
            return
        end
    end
    x = iter(:,k);
end

% Devin Balian 2791430