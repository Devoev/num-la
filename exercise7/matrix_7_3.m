function A = matrix_7_3(n,k)
% Eingabe: n, k >0, k<n
% Ausgabe: nxn Matrix der Bandbreite (k-1)
    kp = [1:k];
    kv = [flip(kp(2:end)), kp];
    k1 = 1./kv;
    on = ones(n,1);
    val = kron(k1,on);
    B = spdiags(val,-(k-1):k-1,n,n);
    D = sparse(diag(1:n));
    A = D*B*D;
    cond(A)
end
