function M = helmholtz_kernel_sp(x, y, k)
if (nargin < 3)
    k = y;
    y = x;
end
r = y-x;
r = sqrt(dot(r,r,2));
M = exp(-1i*k*r)./r;
M(isinf(M)) = 0;
end
