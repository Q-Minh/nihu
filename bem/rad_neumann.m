function [p, resvec] = rad_neumann(model, k, qexc, symm, ng, depth, C)
%RAD_NEUMANN Compute radiation problem with Neumann BC
%   [p, resvec] = rad_neumann(model, ng, depth, k, C, qexc, symm) computes
%   the surface pressure p of a bem radiation problem when the excitation
%   is defined as a velocity constraint.
%  Parameters:
%   model : model mesh
%   k     : wave number
%   qexc  : excitation vector
%   symm  : symmetry parameter, value can be
%              -1 for antysymmetric problems
%               0 for no symmetry (default)
%              +1 for symmetric problems
%           The symmetry plane is the z=0 plane.
%   ng    : number of Gaussian quadrature used for the multipole
%           contribution. The default value is 16 (4x4) for quad elements
%   depth : tree depth, the default value results in a cluster length
%           d = lambda / 4 on the lowest level.
%   C     : Multipole expansion length accuracy parameter, the default
%           value is 3.

% Peter Fiala
% 2009

%% Parameter check
if nargin < 7
    C = 4;
end
if nargin < 6
    D = norm(max(model.Nodes(:,2:4),[],1) - min(model.Nodes(:,2:4),[],1));
    depth = round(log2(k*D/pi));
end
if nargin < 5
    ng = [13 16];
end
if nargin < 4
    symm = 0;
end

%% Gaussian integration points normals and weights
[gcoord, gnorm, w, gind] = geo2gauss(model, ng); % Quadrature config

%% cluster tree
[tree fathersou fatherrec] = clustertree(depth, ecoord, [], symm);  % build cluster tree using elem centers
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
intdata = integpar(tree, k, C, symm);     % determine integration parameters
init_translation(tree, intdata, k, symm); % store translation operators

%% Right hand side
% Excitation in Gaussian points
qw = repmat(qexc.',length(w)/length(qexc),1);
qw = qw(:) .* w;
% multipole contribution
Gqm = mpcont_Gq(ecoord, gcoord, qw, tree, intdata, k, fathergau, fatherrec, symm);
Gq = Gnf * qexc + Gqm;

%% GMRES iteration
restart = 100;
tol = 1e-6;
maxit = 10;
[p,flag,relres,resvec] = mygmres(@(p) gmres_iter_neumann(p, Hnf, ecoord, gcoord, gnorm, w, tree, intdata, k, fathergau, fatherrec, symm),...
    Gq, restart, tol, maxit, M);
resvec = resvec / norm(Gq);
