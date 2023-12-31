function [tree, fathersou, dsrc, fatherrec, drec] = build_regular_clustertree(depth, source, receiver, nearfun)
%BUILD_REGULAR_CLUSTERTREE Build cluster tree

% Peter Fiala
% 2014

% Parameter check
% default arguments

if nargin < 4
    nearfun = @(normdiff)dot(normdiff,normdiff,2) < 4-1e-2;
end

if nargin < 3
    receiver = []; % receiver = source, each node is both source and receiver
end


if isempty(receiver)
    % only source nodes have to be taken into account
    nodes = source;
    % each node is both source and receiver
    is_source = true(size(nodes,1),1);
    is_receiver = is_source;
else
    % sources and receivers are defined separately
    nodes = [source; receiver];
    ns = size(source, 1); % number of sources
    is_source = false(size(nodes,1),1);
    is_source(1:ns) = true;
    is_receiver = ~is_source;
end

% Computing model dimensions
mx = max(nodes,[],1); % largest corner coordinates
r0 = min(nodes,[],1); % smallest corner coordinates
D = max(mx - r0);    % max model size

% allocating structure array
c = cell(depth+1, 1);
tree = struct('diameter', num2cell(D*2.^(0:-1:-depth).'),...
    'coord', c,...
    'father', c,...
    'fatheridx', c,...
    'children', c,...
    'childrenidx', c,...
    'nearfield', c,...
    'interlist', c,...
    'nodsou', c,...
    'nodrec', c,...
    'source', c,...
    'receiver', c);

%
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
    base = r0+d/2; % center of corner cluster
    normnod = round(bsxfun(@minus, nod, base)/d);  % normalized nodes
    % ensure that result is not corrupted by roundoff errors
    normnod = min(max(normnod, 0), 2^l-1);
    [normnod, ~, n] = unique(normnod, 'rows'); % node indices
    tree(iL).coord = bsxfun(@plus, normnod*d, base); % Multipole coordinates
    % Fill nodes, fathers
    A = sparse(n, 1:length(n), 1:length(n));
    A = sort(A,2);
    A = fliplr(full(A(:,any(A,1))));
    if l == depth
        sou = [0; is_source];
        rec = [0; is_receiver];
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
        tree(iL).children = A; % will be used for interaction list of level l+1
        tree(iL+1).father = n;
        sou = [0 tree(iL+1).source.'];
        rec = [0 tree(iL+1).receiver.'];
        tree(iL).receiver = any(reshape(rec(A+1), size(A)), 2);
        tree(iL).source = any(reshape(sou(A+1), size(A)), 2);
    end
end

% fathers
% source fathers
[inod, ~, jnod] = find(tree(end).nodsou);
fathersou(jnod) = inod;
fathersou = fathersou(:);
% receiver fathers
[inod, ~, jnod] = find(tree(end).nodrec);
fatherrec(jnod) = inod;
fatherrec = fatherrec(:);

% nearfield clusters at level 0
l = 0;
iL = l + 1;
tree(iL).nearfield = 1;

%
for l = 1 : depth
    iL = l + 1;
    T = tree(iL);
    Tpar = tree(iL-1);
    
    [sc, ~, rc] = find(Tpar.nearfield);
    [ir, ~, jr] = find(Tpar.children(rc,:));
    [is, ~, js] = find(Tpar.children(sc(ir),:));
    irec = jr(is);
    jsou = js(:);
    
    cond = T.source(jsou) & T.receiver(irec);
    jsou = jsou(cond);
    irec = irec(cond);
    
    normnod = T.coord / T.diameter;
    normdiff = normnod(irec,:) - normnod(jsou,:);
    nearcond = nearfun(normdiff);
    
    jnear = jsou(nearcond);
    inear = irec(nearcond);
    nearfield = sparse(jnear, inear, inear);

    jinter = jsou(~nearcond);
    iinter = irec(~nearcond);
    interlist = sparse(jinter, iinter, iinter);
    
    tree(iL).nearfield = sort(nearfield, 2, 'descend');
    tree(iL).interlist = sort(interlist, 2, 'descend');
end

% compute children, father and near field indices
for l = 1 : length(tree)
    % children indices are not computed on the leaf level
    if l ~= length(tree)
        [m, n] = size(tree(l).children);
        [fath, jsou, chil] = find(tree(l).children);
        trans = (tree(l+1).coord(chil,:) - tree(l).coord(fath,:)) / tree(l+1).diameter;
        childrenidx = trans2idx(trans, 2);
        tree(l).childrenidx = full(sparse(fath, jsou, childrenidx, m, n));
    end
    
    % father indices are not computed on the root level
    if l ~= 1
        trans = (tree(l-1).coord(tree(l).father,:) - tree(l).coord) / tree(l).diameter;
        tree(l).fatheridx = trans2idx(trans, 2);
    end
    
    [sou, jsou, rec] = find(tree(l).interlist);
    trans = (tree(l).coord(rec,:) - tree(l).coord(sou,:)) / tree(l).diameter;
    tree(l).interlistidx = sparse(sou, jsou, trans2idx(trans, 7));
end

dsrc = (source - tree(end).coord(fathersou,:)) / (tree(end).diameter/2);
if ~isempty(receiver)
    drec = (receiver - tree(end).coord(fatherrec,:)) / (tree(end).diameter/2);
end

end
