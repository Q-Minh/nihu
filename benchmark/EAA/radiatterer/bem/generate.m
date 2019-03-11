%%
Le = 5e-2;
[mesh, points] = create_radiatterer(Le);
export_off_mesh(mesh, 'data/radiatterer_05cm_quad.off');
mesh = quad2tria(mesh);
export_off_mesh(mesh, 'data/radiatterer_05cm_tria.off');

fid = fopen('data/points.off', 'w');
fprintf(fid, 'OFF\n');
fprintf(fid, '%d %d 0\n', size(points,1), size(points,1));
fprintf(fid, '%g %g %g\n', points.');
for e = 1 : size(points,1)
    fprintf(fid, '%u %u %u %u\n', [3 e-1 e-1 e-1]);
end
fclose(fid);


fid = fopen('data/points_quad.off', 'w');
fprintf(fid, 'OFF\n');
fprintf(fid, '%d %d 0\n', size(points,1), size(points,1));
fprintf(fid, '%g %g %g\n', points.');
for e = 1 : size(points,1)
    fprintf(fid, '%u %u %u %u %u\n', [4 e-1 e-1 e-1 e-1]);
end
fclose(fid);

