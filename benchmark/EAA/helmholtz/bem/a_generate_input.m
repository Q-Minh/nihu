% a_generate_input.m
% Generate the Input data for the BEM computation

clear;

%% Create data directory as needed
if ~isdir('data')
    mkdir('data');
end

%% Create the mesh and export it
Le = 2e-1;
[mesh, points] = create_helmholtz_bem(Le);
export_off_mesh(mesh, sprintf('data/helmholtz_%02dcm_quad.off', Le*100));
mesh = quad2tria(mesh);
export_off_mesh(mesh, sprintf('data/helmholtz_%02dcm_tria.off', Le*100));

%% Create field points file
fid = fopen('data/points.off', 'w');
fprintf(fid, 'OFF\n');
fprintf(fid, '%d %d 0\n', size(points,1), size(points,1));
fprintf(fid, '%g %g %g\n', points.');
for e = 1 : size(points,1)
    fprintf(fid, '%u %u %u %u\n', [3 e-1 e-1 e-1]);
end
fclose(fid);
