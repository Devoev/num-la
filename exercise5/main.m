% Parse flower.gif
[X,map]=imread("flower.gif");
A=im2double(X,"indexed");

% Calculate SVD and low rank apprximation.
t = 500;
[U,S,V] = svd(A);
At = low_rank(U,S,V,t);

% Show images
imshow(A,map)
figure
imshow(At,map)

% Plot singular values
figure
plot(diag(S))
xlabel('Index $j$', 'Interpreter', 'Latex')
ylabel('Singular value $\sigma_j$', 'Interpreter', 'Latex')