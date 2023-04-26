function [lambda, x] = inv_veciter(A, mu, x0, kmax)
    [n,~] = size(A);
    x = zeros(n,kmax);
    lambda = zeros(1,kmax);

    % Init lambda and x
    x(:,1) = x0;
    lambda(1) = mu;
    y = x0 / norm(x0);
    [L,U,P] = lu(A - mu);

    for k = 2:kmax
        x_k = U\(L\(P*y));
        x(:,k) = x_k;
        y = x_k / norm(x_k);
        lambda(k) = dot(y, A*y);
    end
end