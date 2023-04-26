function [lambda, x] = inv_veciter(A, mu, x0, kmax)
    [n,~] = size(A);
    x = zeros(n,kmax);
    lambda = zeros(1,kmax);

    % Init lambda and x
    x(:,1) = x0 / norm(x0);
    lambda(1) = mu;
    [L,U,P] = lu(A - mu);

    for k = 2:kmax
        x_k = U\(L\(P*x(:,k-1)));
        x_k = x_k / norm(x_k);
        x(:,k) = x_k;
        lambda(k) = dot(x_k, A*x_k);
    end
end