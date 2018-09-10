clear;

%% Create and export geometry
R = 3;
Le = 5*1e-1;
mesh = create_sphere_boundary(R, ceil(R/Le));
mesh = flip_mesh(mesh);
export_off_mesh(mesh, 'data/sphere_quad.off');
mesh_tr = quad2tria(mesh);
export_off_mesh(mesh_tr, 'data/sphere_tria.off');

%% Create and export field points
points = create_sphere_boundary(0.5*R, 4);
points = points.Nodes(:,2:4);

fid = fopen('data/points.off', 'w');
fprintf(fid, 'OFF\n');
fprintf(fid, '%d %d 0\n', size(points,1), size(points,1));
fprintf(fid, '%g %g %g\n', points.');
for e = 1 : size(points,1)
    fprintf(fid, '%u %u %u %u\n', [3 e-1 e-1 e-1]);
end
fclose(fid);

%% Create and export excitations
freqvec = .5 : .5 : 100;
% freqvec = union(freqvec, 52-2 : .2 : 52+2);
% freqvec = union(freqvec, 82.5-2 : .2 : 82.5+2);
% freqvec = union(freqvec, 95.5-2 : .2 : 95.5+2);
nFreqs = length(freqvec);

save data/freqs freqvec

nElem = size(mesh.Elements,1);
exc_const = ones(nElem, nFreqs);
export_excitation(freqvec, exc_const, 'data/quad_const.xct');

nDof = 4 * nElem;
exc_gauss = ones(nDof, nFreqs);
export_excitation(freqvec, exc_gauss, 'data/quad_gauss.xct');
