function mesh = fill_polygon(poly, L)
%FILL_POLYGON  Create 2D TRI-mesh from polygon mesh
%   MESH = FILL_POLYGON(POLY, L)

%%
poly = drop_unused_nodes(poly);
poly = drop_mesh_IDs(poly);

%% Call Mesh2D
hdata.hmax = 1.5*L;
d = poly.Nodes(poly.Elements(:,6),2:3) - poly.Nodes(poly.Elements(:,5),2:3);
q = sqrt(dot(d,d,2));
hdata.edgeh = [(1:length(d))', q];
hdata.fun = @hfun;
hdata.args = {min(q)};
options.dhmax = 1;
options.output = false;
[p,t] = mesh2d(poly.Nodes(:,2:3), poly.Elements(:,5:6), hdata, options);

%% Create NiHu mesh
mesh.Nodes(:,2:3) = p;
mesh.Nodes(:,1) = 1:size(p,1);
mesh.Nodes(:,4) = 0;
mesh.Elements(:,5:7) = t;
mesh.Elements(:,1) = 1:size(t,1);
mesh.Elements(:,2:4) = repmat([ShapeSet.LinearTria.Id 1 1], size(t,1), 1);
[mesh.Materials, mesh.Properties] = default_mat_prop();

end

function h = hfun(x,~,q)
h = q * ones(size(x));
end
