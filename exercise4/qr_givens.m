function [Q, R] = qr_givens(A)
% QR_GIVENS Calculates the QR decomposition of A using givens rotations.

    [n,~] = size(A);
    Q = eye(n);
    R = A;

    for i = 1:n
        for j = 1:i-1
            if R(i,j) == 0
                continue
            end
            G = givens(j, i, R(j,j), R(i,j), n);
            R = G*R;
            Q = Q * G';
        end
    end

end

% Devin Balian 2791430