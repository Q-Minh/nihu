mesh = create_brick_boundary([1 1 1], [10 10 10]);
mesh.Nodes(:,3) = mesh.Nodes(:,3).^2;
mesh.Nodes(:,4) = mesh.Nodes(:,4) + mesh.Nodes(:,3).^2;

mesh = create_sphere_boundary(1, 6);

shapeset = ShapeSet.fromId(mesh.Elements(1,2));

x0 = centnorm(mesh);
p0 = x0(1,:);
it = zeros(size(mesh.Elements,1),1);
for e = 1 : size(mesh.Elements,1)
    coords = mesh.Nodes(mesh.Elements(e,5:8), 2:4).';
    [xi, err, it(e)] = be_invmap(coords, p0, shapeset);
    zeta = xi(end);
    xi = xi(1:end-1);
    [N0, dN0] = shapeset.eval(xi);
    x = N0 * coords.' + cross(dN0(:,:,1) * coords.', dN0(:,:,2) * coords.') * zeta;
    err = norm(x - p0) / norm(p0);
    
    fprintf(1, 'Solution: xi = [%g %f %g], x = [%g %g %g], err: %g\n', xi(1), xi(2), zeta, x(1), x(2), x(3), err);
end

figure;
plot_mesh(mesh, it);
hold on;
plot3(p0(1), p0(2), p0(3), 'r*');
axis equal;
colorbar;