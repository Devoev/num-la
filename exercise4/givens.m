function [G] = givens(p, q, xp, xq, n)
% GIVENS Calculates a givens matrix of size n.
% Inputs:
%   p, q   - Rows and columns effected by the givens matrix.
%   xp, xq - Values of the matrix that should be applied to the givens matrix.
%   n      - Size of the matrix.
% Outputs:
%   G      - The assembled givens matrix.

    G = eye(n);
    r = sqrt(xp^2 + xq^2);
    c = xp/r;
    s = xq/r;
    G(p,p) = c;
    G(q,q) = c;
    G(p,q) = s;
    G(q,p) = -s;
end

% Devin Balian 2791430