function K = helmholtz_3d_hsp_kernel(x, nx, y, ny, k)

rvec = bsxfun(@minus, y, x);
r = sqrt(dot(rvec, rvec, 2));
rdny = dot(rvec, ny, 2) ./ r;
rdnx = (-rvec * nx.') ./ r;
dg = -exp(-1i*k*r) ./ r.^2 .* (1 + 1i*k*r) ./ (4*pi);
ddg = (exp(-1i*k*r).*(- k^2*r.^2 + k*r*2i + 2))./r.^3 / (4*pi);

K = (ddg - dg./r) .* (rdnx .* rdny) - dg./r .* (ny * nx.');

end
