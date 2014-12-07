function M = laplace_matrix_sp(x, y)
if nargin < 2
	y = x;
end
r = y-x;
r = sqrt(dot(r, r, 2));
M = 1./r;
M(isinf(M)) = 0;
end
