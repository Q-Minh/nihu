function g = laplace_3d_slp_kernel(x, ~, y, ~)

rvec = bsxfun(@minus, y, x);
r = sqrt(dot(rvec, rvec, 2));
g = 1 ./ r / (4*pi);

end
