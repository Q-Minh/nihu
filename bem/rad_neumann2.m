function [p, resvec] = rad_neumann2(model, ng, depth, k, C, qexc)
%RAD_NEUMANN2 Compute radiation problem with Neumann BC

%% Gaussian integration points normals and weights
[gcoord, gnorm, w, gind] = geo2gauss(model, ng); % Quadrature config
%% cluster tree
[tree fathersou fatherrec] = clustertree(depth, ecoord);  % build cluster tree using elem centers
print_tree_info(tree);
% Gaussian point fathers
fathergau = zeros(size(gcoord,1),1);
fathergau(gind) = repmat(fathersou,1,size(gind,2));
%% Compute BEM sparse matrices and preconditioner
m = nfij(tree(end).nearfield, tree(end).nodsou, tree(end).nodrec);
[pairs(:,1), pairs(:,2)] = find(m);
[Hnf, Gnf] = bemHG(model, k, 'const', [], pairs);
H = Hnf - 2*pi*speye(size(Hnf));
M = sai(H, tree(end).nodsou);
%% FMBEM integration parameters
intdata = integpar(tree, k, C);     % determine integration parameters
init_translation(tree, intdata, k); % store translation operators
%% Right hand side
% Excitation in Gaussian points
qw = repmat(qexc.',length(w)/length(qexc),1);
qw = qw(:) .* w;
% multipole contribution
Gqm = mpcont_Gq(ecoord, gcoord, qw, tree, intdata, k, fathergau, fatherrec);
Gq = Gnf * qexc + Gqm;
%% GMRES iteration
restart = 100;
tol = 1e-6;
maxit = 10;
[p,flag,relres,resvec] = mygmres(@(p) gmres_iter_neumann(p, Hnf, ecoord, gcoord, gnorm, w, tree, intdata, k, fathergau, fatherrec),...
    Gq, restart, tol, maxit, M);
resvec = resvec / norm(Gq);
