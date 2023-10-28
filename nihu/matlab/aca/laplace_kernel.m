function M = laplace_kernel(x, y)
if nargin < 2
    y = x;
end
dr = zeros(size(x,1), size(y,1));
for d = 1 : size(x,2)
    dr = dr + bsxfun(@minus, x(:,d), y(:,d).').^2;
end
r = sqrt(dr);
M = 1./r;
M(isinf(M)) = 0;
end
