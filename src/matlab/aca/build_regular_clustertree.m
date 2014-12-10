function [tree, fathersou, dsrc, fatherrec, drec] = build_regular_clustertree(depth, source, varargin)
%CLUSTERTREE Build cluster tree 

% Peter Fiala
% 2014

% Parameter check
% default arguments

[~, args] = options(varargin{:});

if numel(args) < 1
    receiver = []; % receiver = source, each node is both source and receiver
else
    receiver = args{1};
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
    nodes = [source; receiver];
    ns = size(source, 1); % number of sources
    issou = false(size(nodes,1),1);
    issou(1:ns) = true;
    isrec = ~issou;
end

% Computing model dimensions
mx = max(nodes,[],1); % largest corner coordinates
r0 = min(nodes,[],1); % smallest corner coordinates
D = max(mx - r0);    % max model size

% allocating structure array
c = cell(depth+1, 1);
tree = struct('diameter', num2cell(D*2.^(0:-1:-depth).'),...
    'coord', c,...
    'father', c, 'fatheridx', c, 'children', c, 'childrenidx', c,...
    'notadmissible', c,...
    'nodsou', c, 'nodrec', c,...
    'source', c, 'receiver', c);

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

% notadmissible clusters at level 0
l = 0;
iL = l + 1;
tree(iL).notadmissible = 1;

%
for l = 1 : depth
    iL = l + 1;
    T = tree(iL);
    T2 = tree(iL-1);
    nC = size(T.coord,1);
    
    so = find(T.source);
    nS = length(so);
    
    normnod = T.coord / T.diameter;
    
    % notadmissible and interaction list matrix
    notadmissible = zeros(nC, 0);
    interlist = zeros(nC, 0);
    for iS = 1 : nS
        iC = so(iS);
        uncle = T2.notadmissible(T.father(iC),:);
        uncle = uncle(uncle ~= 0);
        il = T2.children(uncle,:);
        il = il(il~=0);
        
        % the near field condition
        diff = bsxfun(@minus, normnod(il,:), normnod(iC,:));
        dt = diff(:,3);

        dxy = abs(diff(:,1:2));
        mindt2 = max(dt-1,0).^2;    % minimal dt^2 between two clusters
        maxdt2 = (dt+1).^2;         % maximal dt^2
        mindr2 = dot(max(dxy-1,0), max(dxy-1,0), 2);    % minimal distance
        maxdr2 = dot(dxy+1, dxy+1, 2);  % maximal distance
        
        possible = dt >= 0 & maxdt2 - mindr2 >= 0;
        admissible = mindt2 - maxdr2 > 1;
        
        nadm = il(possible & ~admissible);
        il = il(possible & admissible);
        
        nadm = nadm(T.receiver(nadm));
        il = il(T.receiver(il));
        
        notadmissible(iC, 1:length(nadm)) = nadm;
        interlist(iC, 1:length(il)) = il;
    end
    % clear empty columns
    tree(iL).notadmissible = sparse(notadmissible);
    tree(iL).interlist = sparse(interlist);
end

for l = 1 : length(tree)
    if l ~= length(tree)
        [m, n] = size(tree(l).children);
        [f, j, c] = find(tree(l).children);
        childrenidx = ((tree(l).coord(f,:) - tree(l+1).coord(c,:)) < 0) * [1 2 4]' + 1;
        tree(l).childrenidx = full(sparse(f, j, childrenidx, m, n));
    end
    
    if l ~= 1
        tree(l).fatheridx = ((tree(l-1).coord(tree(l).father,:) - tree(l).coord) < 0) * [1 2 4]' + 1;
    end
end

dsrc = (source - tree(end).coord(fathersou,:)) / (tree(end).diameter/2);
if ~isempty(receiver)
    drec = (receiver - tree(end).coord(fatherrec,:)) / (tree(end).diameter/2);
end


end

