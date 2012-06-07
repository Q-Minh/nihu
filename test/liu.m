function [tree, intdata, p, t] = liu(model, field, k, pinc, qinc, symm)

%% Gaussian integration points normals and weights
ng = [0 16];
[ecoord, gcoord, gnorm, w, gind] = geo2gauss(model, ng); % Quadrature config

%% cluster tree
% determine tree depth
mod = join_models(model, field);
D = norm(max(mod.Nodes(:,2:4),[],1) - min(mod.Nodes(:,2:4),[],1));
depth = round(log2(k*D/pi));
% build tree
ttstart = tic;
points = field.Nodes(:,2:4);
[tree fathersou fatherrec] = clustertree(depth, ecoord, points, symm);  % build cluster tree using elem centers
print_tree_info(tree);
% Gaussian point fathers
fathergau = zeros(size(gcoord,1),1);
fathergau(gind) = repmat(fathersou,1,size(gind,2));
tt = toc(ttstart);

%% FMBEM integration parameters
tmstart = tic;
C = 4;
intdata = integpar(tree, k, C, symm);     % determine integration parameters

init_translation(tree, intdata, k, symm); % store translation operators
tm = toc(tmstart);

%% Right hand side
tistart = tic;
% Excitation in Gaussian points
qw = repmat(qinc.',length(w)/length(qinc),1);
qw = qw(:) .* w;
pw = repmat(pinc.',length(w)/length(pinc),1);
pw = pw(:) .* w;
% multipole contribution
Gqm = mpcont_Gq(points, gcoord, qw, tree, intdata, k, fathergau, fatherrec, symm);
Hpm = mpcont_Hp(points, gcoord, gnorm, pw, tree, intdata, k, fathergau, fatherrec, symm);
p = (Hpm - Gqm)/(4*pi);
ti = toc(tistart);

%% output
t = [tt, tm, ti];
