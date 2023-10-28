function n = node_normals(boundary)
%NODE_NORMALS calculate outward normal vectors in nodes

% Calculate element normals
[~, normals] = centnorm(boundary);

% Number of elements and nodes
nElem = size(boundary.Elements,1);
nNode = size(boundary.Nodes,1);

% Extract elements from model
elements = drop_IDs(boundary);

% Enlarge element matrix
if size(elements,2) < 8
    elements(1,8) = 0;
end

% Extract tria and quad elements
tria = elements(elements(:,2) == 23,[1 (5:7)]);
quad = elements(elements(:,2) == 24,[1 (5:8)]);

% Create sparse matrix
I = [
    reshape(tria(:,2:4).',1,[])
    reshape(quad(:,2:5).',1,[])
    ];
J = [
    reshape(repmat(tria(:,1),1,3).',1,[])
    reshape(repmat(quad(:,1),1,4).',1,[])
    ];

N = sparse(I,J,ones(size(I)),nNode,nElem);

% Calculate node normals
n = N*normals;

% Normalize
n = bsxfun(@times, n, 1./sqrt(dot(n,n,2)));

end

