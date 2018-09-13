function g = helmholtz_3d_slp_singular(x0, ~, corners, k)

[theta, alpha, r] = plane_elem_helper(corners, x0);

g_static = sum(r .* sin(alpha) .* log(tan((alpha + theta)/2) ./ tan(alpha/2)));

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

g_dynamic = w.' * diag(j) *...
    ((exp(-1i*k*r/2) * (-1i*k) .* sinc(k*r/2/pi)));

g =(g_static + g_dynamic) / (4*pi);

end
