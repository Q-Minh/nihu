clear;

%%
L = 5*3*[.3 .4 .5];
Le = 5*1e-1;
mesh = create_brick_boundary(L, ceil(L/Le));
mesh = translate_mesh(mesh, -L/2);
export_off_mesh(mesh, 'data/cube_quad.off');
mesh = quad2tria(mesh);
export_off_mesh(mesh, 'data/cube_tria.off');

points = create_sphere_boundary(1.5*max(L/2), 4);
points = points.Nodes(:,2:4);

fid = fopen('data/points.off', 'w');
fprintf(fid, 'OFF\n');
fprintf(fid, '%d %d 0\n', size(points,1), size(points,1));
fprintf(fid, '%g %g %g\n', points.');
for e = 1 : size(points,1)
    fprintf(fid, '%u %u %u %u\n', [3 e-1 e-1 e-1]);
end
fclose(fid);

