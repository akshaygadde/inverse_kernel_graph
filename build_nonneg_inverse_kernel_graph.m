function W = build_nonneg_inverse_kernel_graph(G)
% Computes a sparse non-negative inverse of Gram matrix G

n = length(G);
nodes = 1:1:n;
W = sparse(n,n);
parfor i = 1:n
    fprintf('node %d\n', i);
    nodes_i = nodes;
    nodes_i(i) = []; 
    G_i = G(nodes_i, nodes_i); % removing i-th row and column from G
    g_i = G(nodes_i,i); % removing the i-th elem from 
    
    options = optimoptions('quadprog','Display','off');
    qpsol = quadprog(G_i, -g_i, [],[],[],[], zeros(n-1,1), [],[], options);
    qpsol(qpsol < 1e-6) = 0;
    temp = zeros(n,1);
    temp(nodes_i) = qpsol;
    W(:,i) = temp;
end

W = max(W,W');


