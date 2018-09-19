function mesh = quad2tria_center(mesh)
%QUAD2TRIA_CENTER  Replace Quad elements with TRIA elements
%   MESH = QUAD2TRIA_CENTER(MESH) replaces each QUAD element in MESH by
%   four TRIA elements. The element IDs are not conserved.
%
% Example:

% last modified: 2018.09.11. FP: first version

narginchk(1, 1);

% keep quad elements
q = mesh.Elements(:,2) == ShapeSet.LinearQuad.Id; % search for quad elements
mesh.Elements = mesh.Elements(q,:);
mesh = drop_unused_nodes(mesh);
mesh = drop_mesh_IDs(mesh);
c = centnorm(mesh);
nN = size(mesh.Nodes,1);
nE = size(mesh.Elements,1);

nodes = [
    mesh.Nodes(:,2:4)
    c
    ];

idx = [1 2; 2 3; 3 4; 4 1];
elems = [];
for e = 1 : 4
    elems = [
        elems
        [mesh.Elements(:,4+idx(e,:)), nN + (1:nE)']
        ];
end

mesh = create_empty_mesh();
mesh.Nodes = [(1:(nN + nE))' nodes];
mesh.Elements(1:4*nE, 1) = 1:4*nE;
mesh.Elements(:, 2) = ShapeSet.LinearTria.Id;
mesh.Elements(:,3:4) = 1;
mesh.Elements(:,5:7) = elems;

end % of function
