function [tree, intdata, p, t] = opti_surf(model, k, pinc, qinc, symm)

%% Gaussian integration points normals and weights
ng = [0 16];
[ecoord, gcoord, gnorm, w, gind] = geo2gauss(model, ng); % Quadrature config

%% cluster tree
% determine tree depth
D = norm(max(model.Nodes(:,2:4),[],1) - min(model.Nodes(:,2:4),[],1));
depth = round(log2(k*D/pi));
% build tree
ttstart = tic;
[tree fathersou fatherrec] = clustertree2(depth, ecoord, [], symm);  % build cluster tree using elem centers
print_tree_info(tree);
% Gaussian point fathers
fathergau = zeros(size(gcoord,1),1);
fathergau(gind) = repmat(fathersou,1,size(gind,2));
tt = toc(ttstart);

%% FMBEM integration parameters
tmstart = tic;
C = 4;
intdata = integpar(tree, k, C, symm);     % determine integration parameters
if nargout > 2
    init_translation(tree, intdata, k, 0); % store translation operators
    tm = toc(tmstart);
    
    %% Right hand side
    tistart = tic;
    % Excitation in Gaussian points
    qw = repmat(qinc.',length(w)/length(qinc),1);
    qw = qw(:) .* w;
    pw = repmat(pinc.',length(w)/length(pinc),1);
    pw = pw(:) .* w;
    % multipole contribution
    Gqm = mpcont_Gq2(ecoord, gcoord, qw, tree, intdata, k, fathergau, fatherrec, symm);
    Hpm = mpcont_Hp2(ecoord, gcoord, gnorm, pw, tree, intdata, k, fathergau, fatherrec, symm);
    p = (Hpm - Gqm)/(4*pi);
    ti = toc(tistart);
    
    %% output
    if nargout > 1
        t = [tt, tm, ti];
    end
end