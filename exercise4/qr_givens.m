function [Q, R] = qr_givens(A, tol, istridiagonal)
% QR_GIVENS Calculates the QR decomposition of A using givens rotations.
% Inputs:
%   A             - The matrix to decompose.
%   tol           - Tolerance to check if elements are zero.
%   istridiagonal - Whether A is of tridiagonal Hessenberg form.
% Outputs:
%   Q             - The orthonorgal matrix.
%   R             - The upper triangular matrix.

    [n,~] = size(A);
    Q = eye(n);
    R = A;

    for i = 2:n
        if istridiagonal
            j_range = i-1;
        else
            j_range = 1:i-1;
        end
        for j = j_range
            if abs(R(i,j)) < tol
                continue
            end
            % Apply givens rotations indirectly
            [Q, R] = givens_apply(Q, R, j, i);
        end
    end
end

% Devin Balian 2791430