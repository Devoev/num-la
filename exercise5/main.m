% Devin Balian 2791430

% Parse flower.gif
[X,map]=imread("flower.gif");
A=im2double(X,"indexed");

% Calculate SVD and low rank apprximation.
t = 300;
[U,S,V] = svd(A);
At = low_rank(U,S,V,t);

% Show images
subplot(1,2,1)
imshow(A,map)
title("Original image")
subplot(1,2,2)
imshow(At,map)
title("Low rank approximation of rank " + t)

% Plot singular values
figure
plot(diag(S))
xlabel('Index $j$', 'Interpreter', 'Latex')
ylabel('Singular value $\sigma_j$', 'Interpreter', 'Latex')