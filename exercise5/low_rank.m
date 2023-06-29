function [At,Ut,St,Vt] = low_rank(U,S,V,t)
% LOW_RANK Calculates a low rank approximation of a given SVD.
% Inputs:
%   U,S,V      - The SVD matrices.
%   t          - The rank of the approximation. Must be lower or equal to the rank of A.
% Outputs:
%   At         - The low rank approximation.
%   Ut,St,Vt   - The low rank SVD matrices.

    [n,m] = size(S);
    if t > min(n, m)
        error("t must be smaller or equal to the rk(A) = " + min(n,m) + "!")
    end
    Ut = U(:,1:t);
    St = S(1:t, 1:t);
    Vt = V(:,1:t);
    At = Ut*St*Vt';
end

% Devin Balian 2791430