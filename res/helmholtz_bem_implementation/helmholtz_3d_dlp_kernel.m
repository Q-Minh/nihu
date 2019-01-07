function g = helmholtz_3d_dlp_kernel(x, ~, y, ny, k)
% HELMHOLTZ_3D_DLP_KERNEL Helmholtz Double layer potential kernel in 3D

rvec = bsxfun(@minus, y, x);
r = sqrt(dot(rvec, rvec, 2));
rdny = dot(rvec, ny, 2) ./ r;
g = -exp(-1i*k*r) ./ r.^2 .* (1 + 1i*k*r) .* rdny / (4*pi);

end
