function g = helmholtz_3d_grad_kernel(x, y, k)
% HELMHOLTZ_3D_DLP_KERNEL Helmholtz Double layer potential kernel in 3D

rvec = bsxfun(@minus, y, x);
r = sqrt(dot(rvec, rvec, 2));
gradr = bsxfun(@times, rvec, 1./ r);
g = bsxfun(@times, -exp(-1i*k*r) ./ r.^2 .* (1 + 1i*k*r) / (4*pi), gradr);

end
