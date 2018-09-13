function g = helmholtz_3d_dlpt_kernel(x, nx, y, ~, k)

rvec = bsxfun(@minus, y, x);
r = sqrt(dot(rvec, rvec, 2));
rdnx = (-rvec * nx.') ./ r;
g = -exp(-1i*k*r) ./ r.^2 .* (1 + 1i*k*r) .* rdnx / (4*pi);

end
