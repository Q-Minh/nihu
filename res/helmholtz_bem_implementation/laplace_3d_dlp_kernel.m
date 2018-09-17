function K = laplace_3d_dlp_kernel(x, ~, y, ny)

rvec = bsxfun(@minus, y, x);
r = sqrt(dot(rvec, rvec, 2));
rdny = dot(rvec, ny, 2) ./ r;
dg = -1 ./ r.^2 ./ (4*pi);

K = dg .* rdny;

end
