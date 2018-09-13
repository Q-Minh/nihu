function g = helmholtz_3d_hsp_singular(x0, ~, corners, k)

[theta, alpha, r] = plane_elem_helper(corners, x0);

IG0 = sum(r .* sin(alpha) .* log(tan((alpha + theta)/2) ./ tan(alpha/2)));
IddG0 = sum((cos(alpha + theta) - cos(alpha)) ./ (r .* sin(alpha)));
g_static = IddG0 + k*k / 2.0 * IG0;

order = 5;
[xi, w] = gaussquad2(order, 3);
[N, dN] = ShapeSet.LinearTria.eval(xi);
y = N * corners;
dyxi = dN(:,:,1) * corners;
dyeta = dN(:,:,2) * corners;
jvec = cross(dyxi, dyeta, 2);
j = sqrt(dot(jvec, jvec, 2));
rvec = bsxfun(@minus, y, x0);
r = sqrt(dot(rvec, rvec, 2));


small = r < 1e-3;
large = ~small;

rl = r(large);
rs = r(small);

k_dynamic = 0 * r;
k_dynamic(large) = (exp(-1i*k*rl).*(1 + 1i*k*rl) - 1 + 1i*k*rl*1i*k.*rl ./ 2.0) ./ rl.^3;
k_dynamic(small) = -1i*k^3 * (...
			1.0 / 3.0 - 1i*k*rs.*(1.0 / 8.0 - 1i*k*rs.*(1.0 / 30.0 - 1i*k*rs.*(1.0 / 144.0 - 1i*k*rs.*(1.0 / 840.0 - 1i*k*rs / 5760.0)))) ...
			);
g_dynamic = w.' * diag(j) * k_dynamic;

g =(g_static + g_dynamic) / (4*pi);

end
