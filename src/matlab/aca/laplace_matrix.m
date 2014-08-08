function M = laplace_matrix(I, J, xr, xs)
if nargin == 3
	xs = xr;
end
x = xr(I,:);
y = xs(J,:);
dx = bsxfun(@minus, x(:,1), y(:,1).');
dy = bsxfun(@minus, x(:,2), y(:,2).');
dz = bsxfun(@minus, x(:,3), y(:,3).');
r = sqrt(dx.^2 + dy.^2 + dz.^2);
M = 1./r;
M(r == 0) = 0;
end
