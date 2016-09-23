close all
clearvars
clc

%%%% m points in R^n
m = 100;
n = 3;

R = orth(randn(n));
t = randn(1,n);
sigma = 0.1; %%% noise level

X = randn(m,n);
Y = X*R + repmat(t,m,1) + sigma*randn(m,n);

[d,Z,transform] = procrustes(Y,X,'scaling',false);

%%%% The output in struct 'transform' includes b (scalar), T (orthogonal
%%%% matrix), and c (matrix representing translation)
%%%% We expect Y = b*X*T + c, so b = 1, T = R, c(1,:) = t
fprintf('scaling: %f\n', transform.b);
fprintf('translation residue: %f\n', norm(t-transform.c(1,:)));
fprintf('orthogonal transform residue: %f\n', norm(R-transform.T));

