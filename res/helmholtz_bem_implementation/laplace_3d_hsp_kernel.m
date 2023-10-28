function K = laplace_3d_hsp_kernel(x, nx, y, ny)

rvec = bsxfun(@minus, y, x);
r = sqrt(dot(rvec, rvec, 2));
rdny = dot(rvec, ny, 2) ./ r;
rdnx = (-rvec * nx.') ./ r;
dg = -1 ./ r.^2 ./ (4*pi);
ddg = 2./r.^3 / (4*pi);

K = (ddg - dg./r) .* (rdnx .* rdny) - dg./r .* (ny * nx.');

end
