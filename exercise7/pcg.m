function [x,k,resvec] = cg(A,b,P,x0,tol,maxit)
% PCG Calculates the solution of the linear system Ax=b using the preconditioned conjugate gradient algorithm (PCG).
% Inputs:
%   A      - System matrix of size (n,n).
%   b      - Right hand side of size(n).
%   P      - Preconditioner matrix.
%   x0     - Initial guess of size(n).
%   tol    - Tolerance.
%   maxit  - Maximum number of iterations.
% Outputs:
%   x      - Approximate solution of size (n).
%   k      - Number of iterations.
%   resvec - Iterated residual errors.

    % Init vars
    resvec = zeros(1,maxit);
    r = b - A*x0;
    x = x0;
    z = P\r;

    for k=1:maxit
        % Calculate new value for rho
        resvec(k) = norm(r);
        if k>1
            rho_old = rho;
        end
        rho = dot(r,z);

        % Determine v
        if k==1
            p = z;
        else
            beta = rho / rho_old;
            p = z + beta*p;
        end

        % Calculate iterate solution x and residual error
        v = A*p;
        alpha = rho / dot(p,v);
        x = x + alpha*p;
        r = r - alpha*v;
        z = P\r;

        % Convergence test
        if resvec(k) < tol
            break
        end
    end
    resvec = resvec(1:k);
end

% Devin Balian 2791430