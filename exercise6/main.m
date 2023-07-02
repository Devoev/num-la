load -ascii 'nos7.mtx'
A=nos7; A=sparse(A(2:end,1),A(2:end,2),A(2:end,3),A(1,1),A(1,2));
load -ascii 'bcsstm10.mtx'
B=nos7; B=sparse(B(2:end,1),B(2:end,2),B(2:end,3),B(1,1),B(1,2));

[n,~] = size(A);
m = 100;
tol = 1e-10;
r0 = ones(n,1);
[V,H] = arnoldi(A,r0,m,tol);

disp("Element h(m+1,m)=" + H(m+1,m))
