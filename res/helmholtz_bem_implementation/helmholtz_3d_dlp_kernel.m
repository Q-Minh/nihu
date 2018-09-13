function g = helmholtz_3d_dlp_kernel(x, ~, y, ny, k)

rvec = bsxfun(@minus, y, x);
r = sqrt(dot(rvec, rvec, 2));
rdny = dot(rvec, ny, 2) ./ r;
g = -exp(-1i*k*r) ./ r.^2 .* (1 + 1i*k*r) .* rdny / (4*pi);

end
