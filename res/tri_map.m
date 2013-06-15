function [x, w] = tri_map(C, N)

[xi, w] = gaussquad2(2*N-1, 3);
[N, dN] = shapefun(xi, 23);

x = N * C;
x_xi = dN(:,:,1) * C;
x_eta = dN(:,:,2) * C;
j = x_xi(:,1).*x_eta(:,2) - x_xi(:,2).*x_eta(:,1);
w = w .* j;

end
