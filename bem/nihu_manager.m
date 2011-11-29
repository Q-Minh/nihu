function nihu_manager(command_file)

load(command_file, 'kvec');

surface_mesh(command_file);

for iK = 1 : numel(kvec)
    surface_sources(command_file, iK);
    surface_tree(command_file, iK);
    surface_sparse(command_file, iK);
    surface_solution(command_file, iK);
end

field_mesh(command_file);

for iK = 1 : numel(kvec)
    field_sources(command_file, iK);
    field_tree(command_file, iK);
    field_sparse(command_file, iK);
    field_solution(command_file, iK);
end

end

%% SURFACE MESH
function surface_mesh(command_file)

tstart = tic;
nihu_log('Opening Surface Mesh...');

load(command_file, 'mesh', 'workdir', 'gauss');
[cent, normal] = centnorm(mesh);
[gs, gn, w, gind] = geo2gauss(mesh, gauss);
save(fullfile(workdir, 'data', 'geometry'),...
    'mesh', 'cent', 'normal', 'gs', 'gn', 'w', 'gind');

nihu_log(sprintf('Ready in %g s.\n', toc(tstart)));
end

%%
function surface_sources(command_file, iK)

tstart = tic;
nihu_log('Computing Surface Source Fields...');

load(command_file, 'workdir', 'src_type', 'symm', 'kvec');
k = kvec(iK);
load(fullfile(workdir, 'data', 'geometry.mat'), 'cent', 'normal');
pis = zeros(size(cent,1),1);    % incident pressure field
qis = zeros(size(cent,1),1);    % incident derivative field
if ~isempty(src_type)
    switch src_type
        case 'point'
            load(command_file, 'r0', 'q0');
            nr = size(r0,1);
            for ir = 1 : nr
                [pisi, qisi] =...
                    incident('point', r0, cent, normal, k, symm);
                pis = pis + q0(ir) * pisi;
                qis = qis + q0(ir) * qisi;
            end
        case 'plane'
            load(command_file, 'dir', 'q0');
            nr = size(dir,1);
            for ir = 1 : nr
                [pisi, qisi] =...
                    incident('plane', dir, cent, normal, k, symm);
                pis = pis + q0(ir) * pisi;
                qis = qis + q0(ir) * qisi;
            end
    end
end
save(fullfile(workdir, 'data', sprintf('incident_k%.3d', iK)),...
    'pis', 'qis');

nihu_log(sprintf('Ready in %g s.\n', toc(tstart)));
end

%%
function surface_tree(command_file, iK)

tstart = tic;
nihu_log('Building Surface Cluster Tree...');

load(command_file, 'workdir', 'kdmin', 'symm', 'kvec');
load(fullfile(workdir, 'data', 'geometry.mat'), 'cent');
k = kvec(iK);
D = max(max(cent,[],1)-min(cent, [], 1));
depth = round(log2(k*D/kdmin));                         % tree depth
[tree, fs, fr] = clustertree(depth, cent, cent, symm);  % build tree
save(fullfile(workdir, 'data', sprintf('surftree_k%.3d', iK)),...
    'tree', 'fs', 'fr');

nihu_log(sprintf('Ready in %g s.\n', toc(tstart)));
end

%%
function surface_sparse(command_file, iK)

tstart = tic;
nihu_log('Assembling Near Field System Matrices...');

load(command_file, 'workdir', 'symm', 'kvec');
k = kvec(iK);
load(fullfile(workdir, 'data', sprintf('surftree_k%.3d', iK)), 'tree');
[i, j] = nfij(tree(end).nearfield, tree(end).nodsou); %#ok<*COLND>
load(fullfile(workdir, 'data', 'geometry'), 'mesh', 'cent');
[Hnf, Gnf] = bemHG(mesh, k, 'const', [], [i j]);
if symm
    [i, j] = nfij(tree(end).imnearfield, tree(end).nodsou);
    if ~isempty(i)
        immesh = reflect_mesh(mesh, [0 0 0], [0 0 1]);
        [imHnf, imGnf] = bemHG(immesh, k, 'const', cent, [i j]);
        Hnf = Hnf + symm * imHnf; %#ok<*NASGU>
        Gnf = Gnf + symm * imGnf;
    end
end
save(fullfile(workdir, 'data', sprintf('surfsparse_k%.3d', iK)),...
    'Hnf', 'Gnf');

nihu_log(sprintf('Ready in %g s.\n', toc(tstart)));
end

%%
function surface_solution(command_file, iK)

tstart = tic;
nihu_log('Computing FMBEM Integration Parameters...');

load(command_file, 'workdir', 'maxit', 'tol', 'symm', 'kvec', 'C');
load(fullfile(workdir, 'data', sprintf('surftree_k%.3d', iK)),...
    'tree', 'fs', 'fr');
load(fullfile(workdir, 'data', 'geometry.mat'),...
    'cent', 'gs', 'gn', 'gind');
k = kvec(iK);
I = integpar(tree, k, C, 1, symm);
[tr, rr, rs, ns] = reltree(tree, cent, fr, gs, fs(gind), gn, symm);
clear tree;
save(fullfile(workdir, 'data', sprintf('surftree_k%.3d', iK)),...
    'tr', 'rr', 'rs', 'ns', 'I', '-append');

nihu_log(sprintf('Ready in %g s.\n', toc(tstart)));

tstart = tic;
nihu_log('Computing Right Hand Side...');

load(command_file, 'qs');
qss = qs;
clear qs;
load(fullfile(workdir, 'data', sprintf('incident_k%.3d', iK)), 'qis');
qss = qss - qis;
clear qis;
load(fullfile(workdir, 'data', sprintf('surfsparse_k%.3d', iK)), 'Gnf');
Gq = Gnf * qss;
clear Gnf
load(fullfile(workdir, 'data', 'geometry'), 'gind', 'w');
Gq = Gq + mpcont_Gq(rr, rs, qss(gind).*w, tr, I, k, fs(gind), fr, symm);

nihu_log(sprintf('Ready in %g s.\n', toc(tstart)));

tstart = tic;
nihu_log('Iterative Solution...');

load(fullfile(workdir, 'data', sprintf('surfsparse_k%.3d', iK)), 'Hnf');
[pss, eps] = my_gmres(...
    @(p) gmres_iter_neumann(p, Hnf, rr, rs, ns, w, gind, tr, I, k, fs, fr, symm),...
    Gq, maxit, tol, speye(size(Hnf)));
clear Hnf I tr gn gs ns rs

load(fullfile(workdir, 'data', sprintf('incident_k%.3d', iK)), 'pis');
pts = pss + pis;
save(fullfile(workdir, 'data', sprintf('solution_k%.3d', iK)),...
    'qss', 'pss', 'pts', 'eps');

nihu_log(sprintf('Ready in %g s.\n', toc(tstart)));
end


%%
function field_mesh(command_file)

tstart = tic;
nihu_log('Opening Field Point Mesh...');

load(command_file, 'workdir', 'field');
points = field.Nodes(:,2:4);
save(fullfile(workdir, 'data', 'geometry'), 'field', 'points', '-append');

nihu_log(sprintf('Ready in %g s.\n', toc(tstart)));
end


%%
function field_sources(command_file, iK)

tstart = tic;
nihu_log('Computing Field Point Source Fields...');

load(command_file, 'workdir', 'kvec', 'src_type', 'symm');
k = kvec(iK);
load(fullfile(workdir, 'data', 'geometry'), 'points');
pif = zeros(size(points,1),1);    % incident pressure field
if ~isempty(src_type)
    switch src_type
        case 'point'
            load(command_file, 'r0', 'q0');
            nr = size(r0,1);
            for ir = 1 : nr
                pif = pif + incident('point', r0, points, [], k, symm);
            end
        case 'plane'
            load(command_file, 'dir', 'q0');
            nr = size(dir,1);
            for ir = 1 : nr
                pif = pif + incident('plane', dir, points, [], k, symm);
            end
    end
end
save(fullfile(workdir, 'data', sprintf('incident_k%.3d', iK)),...
    'pif', '-append');

nihu_log(sprintf('Ready in %g s.\n', toc(tstart)));
end

%%
function field_tree(command_file, iK)

tstart = tic;
nihu_log('Building Field Point Cluster Tree...');

load(command_file, 'workdir', 'kvec', 'symm', 'kdmin');
load(fullfile(workdir, 'data', 'geometry'), 'cent', 'points');
D = max(max([cent; points]) - min([cent; points])); % mesh dimension
depth = round(log2(kvec(iK)*D/kdmin));
[tree, fs, fr] = clustertree(depth, cent, points, symm);
save(fullfile(workdir, 'data', sprintf('fieldtree_k%.3d', iK)),...
    'tree', 'fs', 'fr');

nihu_log(sprintf('Ready in %g s.\n', toc(tstart)));
end

%%
function field_sparse(command_file, iK)

tstart = tic;
nihu_log('Computig Near Field System Matrices...');

load(command_file, 'workdir', 'kvec', 'symm');
k = kvec(iK);
load(fullfile(workdir, 'data', sprintf('fieldtree_k%.3d', iK)), 'tree');
load(fullfile(workdir, 'data', 'geometry'), 'mesh', 'points', 'cent');
[i, j] = nfij(tree(end).nearfield, tree(end).nodsou, tree(end).nodrec);
if ~isempty(i)
    [Hnf, Gnf] = bemHG(mesh, k, 'const', points, [i j]);
else
    Hnf = sparse(size(points,1), size(cent,1));
    Gnf = Hnf;
end
if symm
    [i, j] = nfij(tree(end).imnearfield, tree(end).nodsou, tree(end).nodrec);
    if ~isempty(i)
        immesh = reflect_mesh(mesh, [0 0 0], [0 0 1]);
        [imHnf, imGnf] = bemHG(immesh, k, 'const', points, [i j]);
        Hnf = Hnf + symm * imHnf;
        Gnf = Gnf + symm * imGnf;
    end
end
save(fullfile(workdir, 'data', sprintf('fieldsparse_k%.3d', iK)),...
    'Hnf', 'Gnf');

nihu_log(sprintf('Ready in %g s.\n', toc(tstart)));
end


%%
function field_solution(command_file, iK)

tstart = tic;
nihu_log('Computing FMBEM Integration Parameters...');

load(command_file, 'workdir', 'symm', 'kvec', 'C')
load(fullfile(workdir, 'data', 'geometry'),...
    'points', 'gs', 'gn', 'w', 'gind');
load(fullfile(workdir, 'data', sprintf('fieldtree_k%.3d', iK)),...
    'tree', 'fr', 'fs');
k = kvec(iK);
I = integpar(tree, k, C, 1, symm);
[tr, rr, rs, ns] = reltree(tree, points, fr, gs, fs(gind), gn, symm);
clear tree gn gs

nihu_log(sprintf('Ready in %g s.\n', toc(tstart)));

tstart = tic;
nihu_log('Evaluating FMBEM Integrals...');

load(fullfile(workdir, 'data', sprintf('solution_k%.3d', iK)),...
    'pss', 'qss');
load(fullfile(workdir, 'data', sprintf('fieldsparse_k%.3d', iK)),...
    'Hnf', 'Gnf');
psf = (Hnf * pss - Gnf * qss);
clear Hnf Gnf
psf = psf +...
    mpcont_Hp(rr, rs, ns, pss(gind).*w, tr, I, k, fs(gind), fr, symm);
clear ns pss
psf = psf - mpcont_Gq(rr, rs, qss(gind).*w, tr, I, k, fs(gind), fr, symm);
psf = psf / 4/pi;
clear I tr ns rs rr gind
load(fullfile(workdir, 'data', sprintf('incident_k%.3d', iK)), 'pif');
ptf = pif + psf;    % the total field point pressure
save(fullfile(workdir, 'data', sprintf('solution_k%.3d', iK)),...
    'ptf', 'psf', '-append');

nihu_log(sprintf('Ready in %g s.\n', toc(tstart)));
end

%%
function nihu_log(string)
fprintf(1, string);
end