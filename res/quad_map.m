function [x, w] = quad_map(C, N)

[xi, w] = gaussquad(N);
[xi, eta] = meshgrid(xi, xi);
w = w * w';
xi = [xi(:) eta(:)];
w = w(:);
[N, dN] = shapefun(xi, 24);

x = N * C;
x_xi = dN(:,:,1) * C;
x_eta = dN(:,:,2) * C;
j = x_xi(:,1).*x_eta(:,2) - x_xi(:,2).*x_eta(:,1);
w = w .* j;

end
