function [pf t] = fm_postproc(model, field, k, qexc, pexc, symm, order, depth, C)
%FM_POSTPROC bem post processing with fast multipole algorithm
%   pf = fm_postproc(model, field, k, qexc, pexc, symm, ng, depth, C)
%   [pf t] = fm_postproc(...)

%% Parameter check
if nargin < 9
    C = 4;
end
if nargin < 8
    D = norm(diff(boundingbox(join_meshes(model, field))));
    depth = round(log2(k*D/pi));
end
if nargin < 7
    order = 3;
end
if nargin < 6
    symm = 0;
end

%% Gaussian integration points normals and weights
[gs, gn, w, gind] = geo2gauss(model, order); % Quadrature config

%% cluster tree
ttstart = tic;
points = field.Nodes(:,2:4);
[tree fs fr] = clustertree(depth, gs, points, symm);  % build cluster tree using elem centers
print_tree_info(tree);
% Gaussian point fathers
tt = toc(ttstart);

%% Compute BEM sparse matrices
tnstart = tic;
[i, j] = nfij(tree(end).nearfield, tree(end).nodsou, tree(end).nodrec);
if ~isempty(i)
    [Hnf, Gnf] = bemHG(model, k, 'const', points, [i, j]);
else
    Hnf = sparse([],[],[],size(points,1),length(pexc));
    Gnf = sparse([],[],[],size(points,1),length(pexc));
end
tn = toc(tnstart);

%% FMBEM integration parameters
tmstart = tic;
I = integpar(tree, k, C, symm);     % determine integration parameters
[rt, rr, rs, ns] = reltree(tree, points, fr, gs, fs(gind), gn, symm);
tm = toc(tmstart);

%% Right hand side
tistart = tic;
% multipole contribution
Gqm = mpcont_Gq(rr, rs,     qexc(gind).*w, rt, I, k, fs(gind), fr, symm);
Hpm = mpcont_Hp(rr, rs, ns, pexc(gind).*w, rt, I, k, fs(gind), fr, symm);
pf = Hpm + Hnf * pexc - Gqm - Gnf * qexc;
ti = toc(tistart);

%% output
if nargout > 1
    t = [tt, tn, tm, ti];
end
