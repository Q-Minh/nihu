function g = helmholtz_3d_hsp_nearly_singular(x, nx, corners, k)

% integrate static part analytically
g_stat = laplace_3d_hsp_nearly_singular(x, nx, corners);

% integrate difference Duffyly
[xi, wns] = gaussquad2(20, 4);
[Nns, dNns] = ShapeSet.LinearQuad.eval(xi);

y = Nns * corners([1 2 3 3],:);
ygxi = dNns(:,:,1) * corners([1 2 3 3],:);
ygeta = dNns(:,:,2) * corners([1 2 3 3],:);
jvec = cross(ygxi, ygeta, 2);
jac = sqrt(dot(jvec, jvec, 2));
ny = bsxfun(@times, jvec, 1./jac);
ker = helmholtz_3d_hsp_kernel(x, nx, y, ny, k) ...
    - laplace_3d_hsp_kernel(x, nx, y, ny);
g_dyn = wns.' * diag(jac) * ker;

g = g_stat + g_dyn;

end
