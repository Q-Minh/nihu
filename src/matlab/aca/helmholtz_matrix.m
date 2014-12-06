function M = helmholtz_matrix(I, J, cr, cs, k)
if (nargin < 5)
    k = cs;
    cs = cr;
end
x = cr(I,:);
y = cs(J,:);
dr = bsxfun(@minus, x(:,1), y(:,1)').^2;
dr = dr + bsxfun(@minus, x(:,2), y(:,2)').^2;
dr = dr + bsxfun(@minus, x(:,3), y(:,3)').^2;
dr = sqrt(dr);
M = exp(1i*k*dr)./dr;
M(isinf(M)) = 0;
end
