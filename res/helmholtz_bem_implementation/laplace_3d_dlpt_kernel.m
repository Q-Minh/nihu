function K = laplace_3d_dlpt_kernel(x, nx, y, ~)

rvec = bsxfun(@minus, y, x);
r = sqrt(dot(rvec, rvec, 2));
rdnx = -(rvec * nx.') ./ r;
dg = -1 ./ r.^2 ./ (4*pi);

K = dg .* rdnx;

end
