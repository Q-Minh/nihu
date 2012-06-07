function [tree fathersou fatherrec] = clustertree(depth, source, receiver, symm)
%CLUSTERTREE Build cluster tree of a fast multipole bem model
%   [tree, fs, fr] = clustertree(depth, xs, xr)
%   Builds the cluster tree based on a set of source and receiver nodes.
% Parameters:
%   depth : Depth of the cluster tree (starting with level 0, so the
%           tree consits of depth+1 levels)
%   xs    : Nx3 xyz coordinates of the source nodes
%   xr    : Mx3 xyz coordinates of the receiver nodes. Default is xr = xs
%   symm  : symmetry parameter: -1, 0 (default) or +1
%   tree  : depth+1 long array of tree level structures. Each level L is
%           described by the following fields:
% L.diameter  : Cluster diameter of the current level
% L.coord     : Cx3 xyz cluster centres
% L.father    : Cx1 indices of the father clusters at level l-1
% L.child     : Cx1 indices of the child clusters at level l+1
% L.nearfield : Cx? nearfield cluster indices
% L.interlist : Cx? interaction list cluster indices
% L.nodsou    : Cx? indices of source nodes, only in leaf clusters
% L.nodrec    : Cx? indices of receiver nodes, only in leaf clusters
% L.source    : Cx1 logical, indicates whether a cluster is source
% L.receiver  : Cx1 logical, indicates whether a cluster is receiver
% L.D         : Qx3 xyz unique distances among interaction list entries
% L.Dindex    : Cx? row indices to the matrix D for each interaction list
%               entry
%   fs : Nx1 leaf cluster indices of the source nodes
%   fr : Mx1 leaf cluster indices of the receiver nodes
%
% For symmetric structures (symm = -1 or +1), the tree structures contain
% additional elements:
% L.imnearfield
% L.iminterlist
% L.imDindex

% Peter Fiala
% 2009

%% Parameter check
% default arguments
if nargin < 4
    symm = 0;
end
if nargin < 3
    receiver = []; % receiver = source, each node is both source and receiver
end

if isempty(receiver)
    % only the source nodes have to be taken into account when creating the
    % cluster tree
    nodes = source;
    % each node is both source and receiver
    issou = true(size(nodes,1),1);
    isrec = issou;
else
    % sources and receivers are defined separately
    nodes = [
        source
        receiver
        ];
    ns = size(source, 1); % number of sources
    issou = false(size(nodes,1),1);
    issou(1:ns) = true;
    isrec = ~issou;
end

%% Computing model dimensions
mx = max(nodes,[],1); % largest corner coordinates
mn = min(nodes,[],1); % smallest corner coordinates
if symm
    mn(3) = 0;
end
D = max(mx - mn);    % max model size

%% allocating structure array
c = cell(depth+1,1);
center = mn + D*[1 1 1]/2;
tree = struct('diameter', num2cell(D*2.^(0:-1:-depth).'),...
    'coord', c,...
    'father', c,...
    'child', c,...
    'nearfield', c,...
    'imnearfield', c,...
    'interlist', c,...
    'iminterlist', c,...
    'nodsou', c,...
    'nodrec', c,...
    'source', c,...
    'receiver', c,...
    'D', c,...
    'Dindex', c,...
    'imDindex', c);

%%
r0 = center-D*[1 1 1]/2; % model box lowest corner
% clustering for each level
for l = depth : -1 : 0
    iL = l+1; % level index in the cluster tree structure
    % select nodes to divide into clusters
    if l == depth
        nod = nodes;
    else
        nod = tree(iL+1).coord;
    end
    d = tree(iL).diameter; % current cluster diameter
    base = r0+d/2*[1 1 1]; % center of corner cluster
    normnod = round((nod - repmat(base, size(nod,1),1))/d);  % normalized nodes
    % ensure that result is not corrupted by roundoff errors
    normnod(normnod < 0) = 0;
    normnod(normnod > 2^l-1) = 2^l-1;
    [normnod, ~, n] = unique(normnod, 'rows'); % node indices
    nC = size(normnod,1); % number of clusters on the level
    tree(iL).coord = normnod*d + repmat(base, nC, 1); % Multipole coordinates
    % Fill nodes, fathers
    A = sparse(n, 1:length(n), 1:length(n));
    A = sort(A,2);
    A = fliplr(full(A(:,any(A,1))));
    if l == depth
        sou = [0; issou];
        rec = [0; isrec];
        tree(iL).nodsou = A;
        tree(iL).nodsou(~sou(tree(iL).nodsou+1)) = 0;
        tree(iL).nodrec = A;
        tree(iL).nodrec(~rec(tree(iL).nodrec+1)) = 0;
        if ~isempty(receiver)
            tree(iL).nodrec(logical(tree(iL).nodrec)) = tree(iL).nodrec(logical(tree(iL).nodrec))-ns;
        end
        tree(iL).receiver = any(tree(iL).nodrec,2);
        tree(iL).source = any(tree(iL).nodsou, 2);
    else
        tree(iL).child = A; % will be used for interaction list of level l+1
        tree(iL+1).father = n;
        sou = [0 tree(iL+1).source.'];
        rec = [0 tree(iL+1).receiver.'];
        tree(iL).receiver = any(reshape(rec(A+1), size(A)), 2);
        tree(iL).source = any(reshape(sou(A+1), size(A)), 2);
    end
end

%% fathers
% receiver fathers
[inod, ~, jnod] = find(tree(end).nodrec);
fatherrec(jnod) = inod;
fatherrec = fatherrec(:);
% source fathers
[inod, ~, jnod] = find(tree(end).nodsou);
fathersou(jnod) = inod;
fathersou = fathersou(:);

%% near field and image near field at level 0
l = 0;
iL = l + 1;
tree(iL).nearfield = 1;
if symm
    tree(iL).imnearfield = 1;
end

%%
for l = 1 : depth
    iL = l + 1;
    T = tree(iL);
    T2 = tree(iL-1);
    normnod = T.coord / T.diameter;
    nC = size(tree(iL).coord,1);
    
    so = find(T.source);
    nS = length(so);
    
    % nearfield and interaction list matrix
    nearfield = zeros(nC, 3^3);
    interlist = zeros(nC, 6^3-3^3);
    for iS = 1 : nS
        iC = so(iS);
        progbar(1, nS, iS);
        uncle = T2.nearfield(T.father(iC),:);
        uncle = uncle(uncle ~= 0);
        il = T2.child(uncle,:);
        il = il(il~=0);
        diff = normnod(il,:) - repmat(normnod(iC,:), length(il), 1);
        nf = il(all(round(abs(diff)) <= 1, 2));
        nf = nf(T.receiver(nf));
        il = il(~ismember(il, nf));
        il = il(T.receiver(il));
        nearfield(iC, 1:length(nf)) = nf;
        interlist(iC, 1:length(il)) = il;
    end
    % clear empty columns
    tree(iL).nearfield = nearfield(:,any(nearfield,1));
    tree(iL).interlist = interlist(:,any(interlist,1));
    
    if symm
        % image nearfield and image interaction list matrix
        Trans = diag([1 1 -1]);
        nearfield = zeros(nC, 3^3);
        interlist = zeros(nC, 6^3-3^3);
        for iS = 1 : nS
            iC = so(iS);
            progbar(1, nS, iS);
            uncle = T2.imnearfield(T.father(iC),:);
            uncle = uncle(uncle ~= 0);
            il = T2.child(uncle,:);
            il = il(il~=0);
            diff = normnod(il,:) - repmat(normnod(iC,:)*Trans, length(il), 1);
            nf = il(all(round(abs(diff)) <= 1, 2));
            nf = nf(T.receiver(nf));
            il = setdiff(il, nf);
            il = il(T.receiver(il));
            nearfield(iC, 1:length(nf)) = nf;
            interlist(iC, 1:length(il)) = il;
        end
        % clear empty columns
        tree(iL).imnearfield = nearfield(:,any(nearfield,1));
        tree(iL).iminterlist = interlist(:,any(interlist,1));
    end
end

%% clean the tree
tree = cleantree(tree, symm);

%% distances
tree = clusterdistances(tree, symm);
end

%%
function tree = cleantree(tree, symm)
%cleantree  Clean unneeded levels from the cluster tree
%   tree = cleantree(tree) cuts the levels from the cluster tree where only
%   the structure of a lower level is repeated.

if symm
    lmin = 1;
else
    lmin = 2;
end

ok = true;
while ok && length(tree) > lmin+1
    if size(tree(lmin+2).father,1) == size(tree(lmin+1).coord,1)
        if all(tree(lmin+2).father == (1:size(tree(lmin+1).coord,1)).')
            % copy the 2nd level interaction list to the 3rd level
            tree(lmin+2).interlist = tree(lmin+1).interlist;
            tree(lmin+2).iminterlist = tree(lmin+1).iminterlist;
            % adjust diameter
            tree(lmin+2).diameter = tree(lmin+1).diameter;
            % remove the 2nd level from the tree
            tree = tree([1:lmin, lmin+2:length(tree)]);
        else
            ok = false;
        end
    else
        ok = false;
    end
end
end

%%
function tree = clusterdistances(tree, symm)
%clusterdistances  Compute unique cluster distances in the cluster tree
%   tree = clusterdistances(tree, symm)  computes the unique cluster to
%   cluster distances in all level of the cluster tree.

tree(1).D = zeros(0,3);

depth = length(tree)-1;
for l = 1 : depth
    iL = l+1;
    T = tree(iL);
    nC = length(T.father); % number of clusters in the level
    nIl = size(T.interlist,2); % max size of the interlist
    if symm
        nIl = nIl + max(nIl, size(T.iminterlist,2));
    end
    % allocating space for distance matrices
    D = zeros(nC*nIl, 3); % max size of array of distance vectors
    Dindex = zeros(size(T.interlist)); % distance indices for each cluster pair
    if symm % for symmetric geometries, the distances from the images
        imDindex = zeros(size(T.iminterlist));
        Trans = diag([1 1 -1]); % symmetry operation
    end
    % filling the distance matrices
    ind = 0;
    for iC = 1 : nC % for each cluster
        intlist = T.interlist(iC,:);
        intlist = intlist(intlist~=0);
        nIl = length(intlist);
        if nIl > 0
            ind = max(ind) + (1:nIl);
            % distance between cluster center and interaction clusters
            D(ind,:) = T.coord(intlist,:) - ...
                repmat(T.coord(iC,:), nIl, 1);
            Dindex(iC,1:nIl) = ind;
        end
        if symm
            imintlist = T.iminterlist(iC,:);
            imintlist = imintlist(imintlist~=0);
            nIl = length(imintlist);
            if nIl > 0
                ind = max(ind) + (1:nIl);
                % distance between cluster center and interaction clusters
                D(ind,:) = T.coord(imintlist,:) - ...
                    repmat(T.coord(iC,:) * Trans, nIl, 1);
                imDindex(iC,1:nIl) = ind;
            end
        end
    end
    % cut empty rows from the end
    D = D(1:max(ind),:);
    % compute unique distances
    [D, ~, n] = unique(round(D/T.diameter), 'rows');
    tree(iL).D = D * T.diameter;
    % Dindex should refer to the unique distances
    n2 = [0; n];
    tree(iL).Dindex = n2(Dindex+1);
    if symm % the same operation for the image distances
        tree(iL).imDindex = n2(imDindex+1);
    end
end
end
