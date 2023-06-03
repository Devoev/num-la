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

    if istridiagonal
        for i = 2:n
            j = i-1;
            if abs(R(i,j)) < tol
                continue
            end
            G = givens(j, i, R(j,j), R(i,j), n);
            R = G*R;
            Q = Q*G';
        end
    else
        for i = 1:n
            for j = 1:i-1
                if abs(R(i,j)) < tol
                    continue
                end
                G = givens(j, i, R(j,j), R(i,j), n);
                R = G*R;
                Q = Q*G';
            end
        end
    end
end

% Devin Balian 2791430