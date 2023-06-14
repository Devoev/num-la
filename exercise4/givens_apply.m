function [Q, R] = givens_apply(Q, R, p, q)
% GIVENS_APPLY Applies the givens matrix G to the matrices Q and R without explicitely calculating the matrix products.
% Inputs:
%   Q, R   - The original Q and R matrices.
%   p, q   - Rows and columns effected by the givens matrix.
% Outputs:
%   Q, R   - The modified Q and R matrices.

    % Determine c and s
    [c, s, ~] = givens_coeffs(R(p,p), R(q,p));

    % Modify R and Q
    for k = p:size(R,1)
        Rpk = R(p,k);
        Rqk = R(q,k);
        R(p,k) = c*Rpk + s*Rqk;
        R(q,k) = -s*Rpk + c*Rqk;
    end

    Qp = Q(:,p);
    Qq = Q(:,q);
    Q(:,p) = c*Qp + s*Qq;
    Q(:,q) = -s*Qp + c*Qq;
end

% Devin Balian 2791430