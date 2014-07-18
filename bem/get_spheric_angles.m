function psi = get_spheric_angles(model)
%GET_SPHERIC_ANGLES Outward spheric angles of a surface mesh
%   PSI = GET_SPHERIC_ANGLES(MODEL) Computes the outward spheric angles of
%   the BEM model MODEL.

% Last modified: 15.07.2010.
% Peter Fiala

%% Find model edges and the two neighbor elements for each edge
elements = drop_IDs(model);
edges = get_faces(elements);
% each edge is contained twice, keep unique edges
edg = sort(edges(:,3:4),2);
[edg, m1, n] = unique(edg, 'rows');
nEdg = size(edg,1);
% the matrix el contains the two neighbour elements for each edge
el = zeros(nEdg,2);
el(:,1) = edges(m1,1);
m2 = setdiff(1:length(n), m1);
el(n(m2),2) = edges(m2,1);

%% rearrange matrix el
% For each edge (row), the first element should be directed as the edge
flip = false(nEdg,1);
nvert = mod(elements(:,2),10); % number of vertices for each element
for iEdg = 1 : nEdg
    e1 = el(iEdg,1);
    nv = nvert(e1);
    elem = elements(e1,4+(1:nv));
    g = find(elem == edg(iEdg,1));
    flip(iEdg) = ~(edg(iEdg,2) == elem(mod(g,nv)+1));
end
el(flip,:) = fliplr(el(flip,:));

%% Compute edge angles
[t2, normals] = centnorm(model);
n1 = normals(el(:,1),:);
n2 = normals(el(:,2),:);
d = model.Nodes(edg(:,2),2:4) - model.Nodes(edg(:,1),2:4);
phi = pi + sign(dot(cross(n1, n2, 2), d, 2)) .* acos(dot(n1, n2, 2));

%% Compute spheric angles
nNode = size(model.Nodes,1);
m = sparse(edg(:,1), edg(:,2), phi, nNode, nNode);
ism = m~=0;
psi = sum(m,1)' + sum(m,2) - (sum(ism,1)'+sum(ism,2)-2)*pi;
end
