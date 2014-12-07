function M = helmholtz_matrix(x, y, k)
if (nargin < 3)
    k = y;
    y = x;
end
dr = zeros(size(x,1), size(y,1));
for d = 1 : size(x,2)
    dr = dr + bsxfun(@minus, x(:,d), y(:,d)').^2;
end
dr = sqrt(dr);
M = exp(-1i*k*dr)./dr;
M(isinf(M)) = 0;
end
