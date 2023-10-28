function linmesh = linearise(mesh)
%linmesh convert quadratic tria elements to linear trias

% select quadratic tria elements
qtria = mesh.Elements(:,2) == ShapeSet.QuadraticTria.Id;
qelems = mesh.Elements(qtria, 4+(1:6));
lelems = [
    qelems(:,[1 2 6])
    qelems(:,[2 3 4])
    qelems(:,[4 5 6])
    qelems(:,[2 4 6])
    ];
linmesh = create_empty_mesh();
linmesh.Nodes = mesh.Nodes;
linmesh.Elements(1:size(lelems,1),2) = ShapeSet.LinearTria.Id;
linmesh.Elements(1:size(lelems,1),5:7) = lelems;
linmesh.Elements(1:size(lelems,1),1) = 1:size(lelems,1);
linmesh.Elements(:,3) = 1;
linmesh.Elements(:,4) = 1;
end
