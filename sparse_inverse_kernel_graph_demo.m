clc
clear
close all

%% dependencies
addpath(genpath('data')) % datasets
addpath(genpath('GraphDemos')) % used for simple funcs (e.g. pairwise dist)
addpath(genpath('wgPlot')) % for drawing graphs

%% read data
load('two_moons_1k_0_03.mat'); % data: X (dim x num) and labels: y
y(y==2) = -1; % label for data in class 2 changed to -1 

%% compute graphs
fprintf('Computing the pairwise distance matrix...\n')
D = DistEuclideanPiotrDollar(X',X'); % pairwise Euclidean distances

knn_param = 8;
fprintf('Computing the adj and L of %d-NN graph...\n', knn_param)
knn_mask = sparse(GD_BuildSymmetricKnnGraph(D,knn_param,'dist'));
% sigma in Gaussian kernel = mean(k-th nearest neighbor distance)
sigma = 10*mean(max(D.*knn_mask, [], 2));
fprintf('RBF kernel similarity with sigma = %.4f\n', sigma)
G = exp(-D/(2*sigma^2)); % Gram matrix
W_knn = G .* knn_mask;
L_knn = compute_laplacians(W_knn);

fprintf('Computing the adj and L nonneg inverse kernel...\n')
W_nnik = build_nonneg_inverse_kernel_graph(G);
L_nnik = compute_laplacians(W_nnik);

%% comparing two graphs
num_eigs = 4;
fprintf('Computing the first %d eigs...\n', num_eigs)
[V_knn,lambda_knn] = sorted_eigs(L_knn,num_eigs);
[V_nnik,lambda_nnik] = sorted_eigs(L_nnik,num_eigs);

fprintf('Clustering with 2nd eigvec...\n')
fprintf('with knn graph <v,y>=%.4f\n', V_knn(:,2)'*y)
fprintf('with nnik graph <v,y>=%.4f\n', V_nnik(:,2)'*y)

%% plots


