function [V,lambda] = sorted_eigs(L,m,mode)
% computes the smallest m eigenvals and eigenvecs arranged in ascending order

[V,lambda] = eigs(L,m,'sm');
[lambda, idx] = sort(diag(lambda),'ascend');
V = V(:, idx);