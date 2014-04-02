function M = laplace_matrix(I, J, c)
x = c(I,:);
y = c(J,:);
dx = bsxfun(@minus, x(:,1), y(:,1)');
dy = bsxfun(@minus, x(:,2), y(:,2)');
dz = bsxfun(@minus, x(:,3), y(:,3)');
r = sqrt(dx.^2 + dy.^2 + dz.^2);
M = 1./r;
M(isinf(M)) = 0;
end
