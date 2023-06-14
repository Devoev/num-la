function [G] = givens(p, q, xp, xq, n)
% GIVENS Calculates a givens matrix of size n.
% Inputs:
%   p, q   - Rows and columns effected by the givens matrix.
%   xp, xq - Values of the matrix that should be applied to the givens matrix.
%   n      - Size of the matrix.
% Outputs:
%   G      - The assembled (sparse) givens matrix.

    % Determine r, c and s
    [c, s, r] = givens_coeffs(xp, xq);

    % Assamble sparse matrix
    i = [p, p, q, q];
    j = [p, q, p, q];
    v = [c - 1, s, -s, c - 1]; % Use c-1, because the c's are on the diagonal
    G = speye(n) + sparse(i, j, v, n, n);
end

% Devin Balian 2791430