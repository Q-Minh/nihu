clear;

%% Export radiatterer meshes
Le = 10e-2/3;
[qmesh, points] = create_radiatterer(Le);
export_off_mesh(qmesh, sprintf('data/radiatterer_%03dmm_quad.off', floor(1000*Le)));
tmesh = quad2tria(qmesh);
export_off_mesh(tmesh, sprintf('data/radiatterer_%03dmm_tria.off', floor(1000*Le)));

%% Export field point meshes
% print field points with quad semantics
fid = fopen('data/radiatterer_points_tria.off', 'w');
fprintf(fid, 'OFF\n');
fprintf(fid, '%d %d 0\n', size(points,1), size(points,1));
fprintf(fid, '%g %g %g\n', points.');
for e = 1 : size(points,1)
    fprintf(fid, '%u %u %u %u\n', [3 e-1 e-1 e-1]);
end
fclose(fid);

%% export field plane
% print field points with tria semantics
fid = fopen('data/radiatterer_points_quad.off', 'w');
fprintf(fid, 'OFF\n');
fprintf(fid, '%d %d 0\n', size(points,1), size(points,1));
fprintf(fid, '%g %g %g\n', points.');
for e = 1 : size(points,1)
    fprintf(fid, '%u %u %u %u %u\n', [4 e-1 e-1 e-1 e-1]);
end
fclose(fid);

field = translate_mesh(create_slab([4 4], [4 4]/Le), [-.75, -1 .5]);
hold on;
plot_mesh(field);
x = centnorm(field);

bricks = {
    [0.2, 1.4], [0.0, 1.2], [0.4, 1.2]
    [2.3, 2.5], [0.2, 0.4], [0.5, 0.7]
    [2.3, 2.5], [0.7, 0.8], [1.0, 1.1]
    [1.6, 2.3], [0.1, 1.2], [0.2, 1.5]
    [0.2, 1.0], [1.4, 1.8], [0.0, 1.7]
    [1.2, 1.4], [1.4, 2.0], [1.6, 1.7]
    [1.6, 1.8], [1.4, 2.0], [1.6, 1.7]
    [2.0, 2.2], [1.4, 2.0], [1.6, 1.7]
    [2.4, 2.5], [1.4, 2.0], [1.6, 1.7]
    [1.2, 1.4], [1.4, 2.0], [1.2, 1.4]
    [1.6, 1.8], [1.4, 2.0], [1.2, 1.4]
    [2.0, 2.2], [1.4, 2.0], [1.2, 1.4]
    [2.4, 2.5], [1.4, 2.0], [1.2, 1.4]
    [1.2, 1.4], [1.4, 2.0], [0.6, 1.0]
    [1.6, 1.8], [1.4, 2.0], [0.6, 1.0]
    [2.0, 2.2], [1.4, 2.0], [0.6, 1.0]
    [2.4, 2.5], [1.4, 2.0], [0.6, 1.0]
    [1.2, 1.4], [1.4, 2.0], [0.0, 0.4]
    [1.6, 1.8], [1.4, 2.0], [0.0, 0.4]
    [2.0, 2.2], [1.4, 2.0], [0.0, 0.4]
    [2.4, 2.5], [1.4, 2.0], [0.0, 0.4]
    };

L = [2.5 2.0 1.7];
outside = any(x < 0, 2) | x(:,1) > L(1) | x(:,2) > L(2) | x(:,3) > L(3);
for ib = 1 : length(bricks)
    outside = outside |...
        (x(:,1) > bricks{ib,1}(1) & x(:,1) < bricks{ib,1}(2) & ...
        x(:,2) > bricks{ib,2}(1) & x(:,2) < bricks{ib,2}(2) & ...
        x(:,3) > bricks{ib,3}(1) & x(:,3) < bricks{ib,3}(2));
end

export_off_mesh(field, sprintf('data/radi_plane_%03dmm_quad.off', Le*1000));
