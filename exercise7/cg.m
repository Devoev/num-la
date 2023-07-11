function [x,k,resvec] = cg(A,b,x0,tol,maxit)
% CG Calculates the solution of the linear system Ax=b using the conjugate gradient algorithm (CG).
% Inputs:
%   A      - System matrix of size (n,n).
%   b      - Right hand side of size(n).
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

    for k=1:maxit
        % Calculate new value for rho
        resvec(k) = norm(r);
        rho = resvec(k)^2;

        % Determine v
        if k==1
            p = r;
        else
            beta = rho / resvec(k-1)^2;
            p = r + beta*p;
        end

        % Calculate iterate solution x and residual error
        v = A*p;
        alpha = rho / dot(p,v);
        x = x + alpha*p;
        r = r - alpha*v;

        % Convergence test
        if rho < tol
            break
        end
    end
    resvec = resvec(1:k);
end

% Devin Balian 2791430