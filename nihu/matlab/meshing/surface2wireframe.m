function wireframe = surface2wireframe(mesh, angletol)
%SURFACE2WIREFRAME 
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Extract surface TRIA and QUAD elements
% mesh.Elements = mesh.Elements(floor(mesh.Elements(:,2)/10) == 2,:);

wireframe = get_boundary(mesh);

% Compute element normals
[~, normal] = centnorm(mesh);

% Extract model edges
mesh.Elements = drop_IDs(mesh);
edges = get_faces(mesh.Elements);
% Each edge is contained twice, keep unique edges
edg = sort(edges(:,3:4),2);
[edg, m1, n] = unique(edg, 'rows');
nEdg = size(edg,1);
% the matrix el contains the two neighbour elements for each edge
el = zeros(nEdg,2);
el(:,1) = edges(m1,1);
m2 = setdiff(1:length(n), m1);
el(n(m2),2) = edges(m2,1);

% angle compute
bind = find(el(:,2) ~= 0);
el1 = el(bind,:);
angles = 1-dot(normal(el1(:,1),:), normal(el1(:,2),:), 2);
wf.Elements(:,5:6) = edg(bind(angles > 1-cos(angletol)),:);
wf.Nodes = mesh.Nodes;
nEl = size(wf.Elements,1);
wf.Elements(:,1) = 1:size(wf.Elements,1);
wf.Elements(:,2:4) = repmat([ShapeSet.LinearLine.Id 1 1], nEl,1);
[wf.Materials, wf.Properties] = default_mat_prop();

%
wireframe = join_meshes(wireframe, wf);
end
