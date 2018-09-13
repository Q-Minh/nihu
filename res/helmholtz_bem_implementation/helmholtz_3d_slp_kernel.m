function g = helmholtz_3d_slp_kernel(x, ~, y, ~, k)

rvec = bsxfun(@minus, y, x);
r = sqrt(dot(rvec, rvec, 2));
g = exp(-1i*k*r) ./ r / (4*pi);

end
