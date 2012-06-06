function [pf t] = fm_postproc(model, field, k, qexc, pexc, symm, order, depth, C)
%FM_POSTPROC bem post processing with fast multipole algorithm
%   pf = fm_postproc(model, field, k, qexc, pexc, symm, ng, depth, C)
%   [pf t] = fm_postproc(...)

%% Parameter check
if nargin < 9
    C = 4;
end
if nargin < 8
    mod = join_meshes(model, field);
    D = norm(max(mod.Nodes(:,2:4),[],1) - min(mod.Nodes(:,2:4),[],1));
    depth = round(log2(k*D/pi));
end
if nargin < 7
    order = 3;
end
if nargin < 6
    symm = 0;
end

%% Gaussian integration points normals and weights
[gcoord, gnorm, w, gind] = geo2gauss(model, order); % Quadrature config

%% cluster tree
ttstart = tic;
points = field.Nodes(:,2:4);
[tree fathersou fatherrec] = clustertree(depth, gcoord, points, symm);  % build cluster tree using elem centers
print_tree_info(tree);
% Gaussian point fathers
fathergau = zeros(size(gcoord,1),1);
fathergau(gind) = repmat(fathersou,1,size(gind,2));
tt = toc(ttstart);

%% Compute BEM sparse matrices
tnstart = tic;
[i, j] = nfij(tree(end).nearfield, tree(end).nodsou, tree(end).nodrec);
if ~isempty(i)
    [Hnf, Gnf] = bemHG(model, k, 'const', points, [i, j]);
else
    Hnf = 0;
    Gnf = 0;
end
tn = toc(tnstart);

%% FMBEM integration parameters
tmstart = tic;
intdata = integpar(tree, k, C, symm);     % determine integration parameters
init_translation(tree, intdata, k, symm); % store translation operators
tm = toc(tmstart);

%% Right hand side
tistart = tic;
% Excitation in Gaussian points
qw = repmat(qexc.',length(w)/length(qexc),1);
qw = qw(:) .* w;
pw = repmat(pexc.',length(w)/length(pexc),1);
pw = pw(:) .* w;
% multipole contribution
Gqm = mpcont_Gq(points, gcoord, qw, tree, intdata, k, fathergau, fatherrec, symm);
Hpm = mpcont_Hp(points, gcoord, gnorm, pw, tree, intdata, k, fathergau, fatherrec, symm);
pf = Hpm + Hnf * pexc - Gqm - Gnf * qexc;
ti = toc(tistart);

%% output
if nargout > 1
    t = [tt, tn, tm, ti];
end
