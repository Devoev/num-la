load -ascii 'nos7.mtx'
A=nos7; A=sparse(A(2:end,1),A(2:end,2),A(2:end,3),A(1,1),A(1,2));
load -ascii 'bcsstm10.mtx'
B=bcsstm10; B=sparse(B(2:end,1),B(2:end,2),B(2:end,3),B(1,1),B(1,2));

%% (a) Arnoldi test
[n,~] = size(A);
m = 100;
tol = 1e-8;
r0 = ones(n,1);
[V,H] = arnoldi(A,r0,m,tol);

disp("Element h(m+1,m)=" + H(m+1,m))
Vm = V(:,1:m);
Hm = H(1:m,:);
err = norm(A*Vm - V*H);
assert(err < tol, "Arnoldi method doesn't converge for given matrix A and given tolerance=" + tol + "! Error=" + err)
% TODO: Fix 2nd error test
%err = norm(Vm'*A*Vm - Hm);
%err = norm(Vm'*Vm - speye(m,m));
%assert(err < tol, "Arnoldi method doesn't converge for given matrix A and given tolerance=" + tol + "! Error=" + err)

%% (b) FOM test
[n,~] = size(B);
tol = 1e-9;
kmax = 500;
b = ones(n,1);
[x,iter] = fom(B,b,kmax,tol);

err = norm(B*x - b)/norm(b);
assert(err < tol, "FOM doesn't converge for given matrix B and given tolerance=" + tol + "! Error=" + err)

x_exact = B\b;
err_k = vecnorm(iter - x_exact) / norm(x_exact);
semilogy(err_k)
title("FOM convergence study")
xlabel('Iteration number $k$', 'Interpreter', 'Latex')
ylabel('Relative $L^2$ error $\frac{||x^{(k)} - x||_2}{||x||_2}$', 'Interpreter', 'Latex')
yline(tol, "r--")
legend("Relative error", "Error tolerance")
grid on
