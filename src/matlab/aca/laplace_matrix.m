function M = laplace_matrix(x, y)
if nargin < 2
	y = x;
end
dr = bsxfun(@minus, x(:,1), y(:,1).').^2;
dr = dr + bsxfun(@minus, x(:,2), y(:,2).').^2;
dr = dr + bsxfun(@minus, x(:,3), y(:,3).').^2;
r = sqrt(dr);
M = 1./r;
M(isinf(M)) = 0;
end
