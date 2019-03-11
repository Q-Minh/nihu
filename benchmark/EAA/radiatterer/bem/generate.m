clear;
%%
Le = 5e-2;
[qmesh, points] = create_radiatterer(Le);
export_off_mesh(qmesh, 'data/radiatterer_05cm_quad.off');
tmesh = quad2tria(qmesh);
export_off_mesh(tmesh, 'data/radiatterer_05cm_tria.off');

fid = fopen('data/radiatterer_points_tria.off', 'w');
fprintf(fid, 'OFF\n');
fprintf(fid, '%d %d 0\n', size(points,1), size(points,1));
fprintf(fid, '%g %g %g\n', points.');
for e = 1 : size(points,1)
    fprintf(fid, '%u %u %u %u\n', [3 e-1 e-1 e-1]);
end
fclose(fid);


fid = fopen('data/radiatterer_points_quad.off', 'w');
fprintf(fid, 'OFF\n');
fprintf(fid, '%d %d 0\n', size(points,1), size(points,1));
fprintf(fid, '%g %g %g\n', points.');
for e = 1 : size(points,1)
    fprintf(fid, '%u %u %u %u %u\n', [4 e-1 e-1 e-1 e-1]);
end
fclose(fid);

