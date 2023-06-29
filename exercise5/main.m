% Devin Balian 2791430

% Parse flower.gif
[X,map]=imread("flower.gif");
A=im2double(X,"indexed");

% Calculate SVD and low rank apprximation.
t = 50;
[U,S,V] = svd(A);
[At,Ut,St,Vt] = low_rank(U,S,V,t);

% Show images
subplot(1,2,1)
imshow(A,map)
title("Original image")
subplot(1,2,2)
imshow(At,map)
title("Low rank approximation of rank " + t)

% Plot singular values
figure
semilogy(diag(S))
xline(t, "r--")
xlabel('Index $j$', 'Interpreter', 'Latex')
ylabel('Singular value $\sigma_j$', 'Interpreter', 'Latex')
legend("Singular values", "Low rank cut-off at " + t)

% Calculate storage (only nonzero elements must be stored when using sparse matrices)
storage_full = nnz(S) + nnz(U) + nnz(V);
storage_low_rank = nnz(St) + nnz(Ut) + nnz(Vt);
storage_saved_rel = (storage_full - storage_low_rank) / storage_full;
disp("Saved storage by low rank approximation: " + storage_saved_rel + "%")
