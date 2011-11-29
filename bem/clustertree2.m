function [tree fathersou fatherrec] = clustertree2(depth, source, receiver, symm)
%CLUSTERTREE2 Build cluster tree of a fast multipole bem model
%   [TREE, FATHERSOU, FATHERREC] = CLUSTERTREE(DEPTH, SOURCE, RECEIVER)
%   Builds the cluster tree based on a set of source and receiver nodes.
% Parameters:
%   DEPTH    : Depth of the cluster tree (starting with level 0, so the
%              tree consits of DEPTH+1 levels)
%   SOURCE   : N x 3 matrix containiong the xyz coordinates of the source
%              nodes
%   RECEIVER : M x 3 matrix containing the xyz coordinates of the receiver
%              nodes. Optional, if omitted RECEIVER = SOURCE.
%   TREE     : DEPTH+1 long array of tree level structures. Each level L is
%              described by the following fields:
% L.diameter  : Cluster diameter of the current level
% L.coord     : C x 3 matrix, containing xyz coordinates of cluster centres
%               (multipoles) at the current level
% L. father   : C x 1 matrix containing the indices of the father clusters
%               at a higher level
% L.nearfield : C x ? matrix contatining the nearfield cluster indices
% L.interlist : C x ? matrix containing the interaction list cluster
%               indices
% L.nodsou    : Used only at leaf level. C x ? matrix containing indices of
%               SOURCE nodes for each cluster
% L.nodrec    : Used only at leaf level. C x ? matrix containing indices of
%               RECEIVER nodes for each cluster
% L.source    : Cx1 vector, indicates whether a cluster is source or not
% L.receiver  : Cx1 vector, indicates whether a cluster is receiver or not
% D           : Q x 3 vector: unique distances between interaction list
%               elements
% Dindex      : C x ? matrix, containing row indices to the matrix D for
%               each cluster pair in the interaction lists
%   FATHERSOU : N x 1 matix containing the cluster indices of the source
%               nodes (indexing rows of L(end).Coord)
%   FATHERREC : M x 1 matix containing the cluster indices of the receiver
%               nodes (indexing rows of L(end).Coord)

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
r0 = center-D*[1 1 1]/2; % extended model box lowest corner
% clustering for each level
for l = depth : -1 : 0
    iL = l+1; % level index in the cluster tree
    % select nodes to divide into clusters
    if l == depth
        nod = nodes;
    else
        nod = tree(iL+1).coord;
    end
    d = tree(iL).diameter; % current cluster diameter
    base = r0+d/2*[1 1 1]; % center of corner cluster
    normnod = (nod - repmat(base, size(nod,1),1))/d;  % normalized nodes
    normnod = round(normnod);
    normnod(normnod < 0) = 0;
    normnod(normnod > 2^l-1) = 2^l-1;
    [normnod, m, n] = unique(round(normnod), 'rows'); % node indices
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


%% near field and image near field at level 0
l = 0;
iL = l + 1;
tree(iL).nearfield = 1;
tree(iL).imnearfield = 1;

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
    interlist = zeros(nC, 7^3-3^3);
    for iS = 1 : nS
        iC = so(iS);
        progbar(1, nS, iS);
        % if the given cluster is not source then continue
        uncle = T2.nearfield(T.father(iC),:);
        uncle = uncle(uncle ~= 0);
        il = T2.child(uncle,:);
        il = il(il~=0);
        diff = normnod(il,:) - repmat(normnod(iC,:), length(il), 1);
        nf = il(all(round(abs(diff)) <= 1, 2));
        nf = nf(T.receiver(nf));
        il = setdiff(il, nf);
        il = il(T.receiver(il));
        nearfield(iC, 1:length(nf)) = nf;
        interlist(iC, 1:length(il)) = il;
    end
    % clear empty columns
    tree(iL).nearfield = nearfield(:,any(nearfield,1));
    tree(iL).interlist = interlist(:,any(interlist,1));
    
    % image nearfield and image interaction list matrix
    Trans = diag([1 1 -1]);
    nearfield = zeros(nC, 3^3);
    interlist = zeros(nC, 7^3-3^3);
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

%% fathers
% receiver fathers
[inod, thrash, jnod] = find(tree(end).nodrec);
fatherrec(jnod) = inod;
% source fathers
[inod, thrash, jnod] = find(tree(end).nodsou);
fathersou(jnod) = inod;

%% clean the tree
tree = cleantree(tree);

%% distances
tree = clusterdistances(tree, symm);
end

function tree = cleantree(tree)
%cleantree  Clean unneeded levels from the cluster tree
%   tree = cleantree(tree) cuts the levels from the cluster tree where only
%   the structure of a lower level is repeated.

ok = true;
while ok && length(tree) > 3
    if size(tree(3+1).father,1) == size(tree(2+1).coord,1)
        if ~any(tree(3+1).father - (1:size(tree(2+1).coord,1)).')
            % copy the 2nd level interaction list to the 3rd level
            tree(3+1).interlist = tree(2+1).interlist;
            tree(3+1).iminterlist = tree(2+1).iminterlist;
            % remove the 2nd level from the tree
            ind = [1:2, 4:length(tree)];
            tree = tree(ind);
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
%clusterdistances  Compute uniqu cluster distances in the cluster tree
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
    [D, thrash, n] = unique(round(D/T.diameter), 'rows');
    tree(iL).D = D * T.diameter;
    % Dindex should refer to the unique distances
    n2 = [0; n];
    tree(iL).Dindex = n2(Dindex+1);
    if symm % the same operation for the image distances
        tree(iL).imDindex = n2(imDindex+1);
    end
end
end
