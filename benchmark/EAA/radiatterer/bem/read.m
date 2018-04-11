clear;

%%
Le = 1e-1;
[mesh, points] = create_radiatterer(Le);
export_off_mesh(mesh, 'data/radiatterer_10cm_quad.off');
mesh = quad2tria(mesh);
export_off_mesh(mesh, 'data/radiatterer_10cm_tria.off');

fid = fopen('data/points.off', 'w');
fprintf(fid, 'OFF\n');
fprintf(fid, '%d %d 0\n', size(points,1), size(points,1));
fprintf(fid, '%g %g %g\n', points.');
for e = 1 : size(points,1)
    fprintf(fid, '%u %u %u %u\n', [3 e-1 e-1 e-1]);
end
fclose(fid);

%%
c = 340;
fmax = min(mesh_kmax(mesh)) * c / 2/pi;

%%
res = importdata('data/response.res');
k = res(1);
ps = complex(res(2:2:end), res(3:2:end));
ps = reshape(ps, N, []).';
figure;
plot_mesh(mesh, mean(abs(ps), 2));

