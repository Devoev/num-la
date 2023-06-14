function [c, s, r] = givens_coeffs(xp, xq)
% GIVENS_COEFFS Calculates the c, s and r coefficients for the givens matrix.
% Inputs:
%   xp, xq  - Values of the matrix that should be applied to the givens matrix.
% Outputs:
%   c, s, r - Calculated coefficients.

    r = sqrt(xp^2 + xq^2);
    c = xp/r;
    s = xq/r;
end

% Devin Balian 2791430