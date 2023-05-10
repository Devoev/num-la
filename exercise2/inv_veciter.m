function [lambda, x] = inv_veciter(A, mu, x0, kmax)
    [n,~] = size(A);
    x = zeros(n,kmax);
    lambda = zeros(1,kmax);

    % Init lambda and x
    y = x0 / norm(x0);
    [L,U,P] = lu(A - mu*eye(n));

    for k = 1:kmax
        x_k = U\(L\(P*y));
        y = x_k / norm(x_k);
        x(:,k) = y;
        lambda(k) = dot(y, A*y);
    end
end