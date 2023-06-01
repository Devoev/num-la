function [G] = givens(A, p, q, j)
% GIVENS Calculates the givens matrix for A that sets the element A(q,j) to zero.

    [n, ~] = size(A);
    G = eye(n);
    xp = A(p,j);
    xq = A(q,j);
    r = sqrt(xp^2 + xq^2);
    c = xp/r;
    s = xq/r;
    G(p,p) = c;
    G(q,q) = c;
    G(p,q) = s;
    G(q,p) = -s;
end

% Devin Balian 2791430