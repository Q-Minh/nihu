function K = laplace_3d_grad_kernel(x, ~, y, ~)

rvec = bsxfun(@minus, y, x);
r = sqrt(dot(rvec, rvec, 2));
gradr = bsxfun(@times, rvec, 1./r);
dg = -1 ./ r.^2 ./ (4*pi);

K = bsxfun(@times, dg, gradr);

end
