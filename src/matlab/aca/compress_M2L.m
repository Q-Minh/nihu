function [U, V, r] = compress_M2L(M2L, eps)
%COMPRESS_M2L Compress M2L matrices with ACA

N = size(M2L,1);
K = size(M2L,3);
U = zeros(N,0,K);
V = zeros(N,0,K);
r = zeros(K,1);
for k = 1 : K
    [u, v] = lowrank_approx(M2L(:,:,k), [N, N], eps, N);
    r(k) = size(u,2);
    U(:,1:r(k),k) = u;
    V(:,1:r(k),k) = v;
end

end
