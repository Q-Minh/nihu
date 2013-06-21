function [x, w] = quadrature_map(xi, w, C, type)

[N, dN] = shapefun(xi, type);

x = N * C;
x_xi = dN(:,:,1) * C;
x_eta = dN(:,:,2) * C;
j = x_xi(:,1).*x_eta(:,2) - x_xi(:,2).*x_eta(:,1);
w = w .* j;

end
