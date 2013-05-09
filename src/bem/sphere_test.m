%% compile C++ mex code
clc;
mex -v CXXFLAGS="\$CXXFLAGS -std=c++11 -O3" Boonen13.mex.cpp -I../../../eigen -output Boonen13

%% Build a good little mesh
R = 1;
mesh = create_sphere_boundary(1, 4);
[~, ~, w] = geo2gauss(mesh, 0);
k = min(mesh_kmax(mesh)) / 2;

%% Collocational matrices C++ version
qmesh = quadratise(mesh);
qmesh.Nodes(:,2:4) = bsxfun(@times, qmesh.Nodes(:,2:4), 1./sqrt(dot(qmesh.Nodes(:,2:4), qmesh.Nodes(:,2:4), 2)));
tic;
[nodes, elements] = extract_Boonen_mesh(qmesh);
[G, H] = Boonen13(nodes, elements, k);
tBooni = toc;

%%
tic;
[nodes, elements] = extract_bem_mesh(mesh);
[HH, GG] = bemHG(mesh, k, 'lin');
tOldSchool = toc;

%% Solve
M = ibemRHS(qmesh);
q = ones(size(G,1),1);
p = (H - .5*M) \ (G * q);
q = ones(size(GG,1),1);
pp = (HH - .5*eye(size(GG,1))) \ (GG * q);

%% Analytical
p_anal = -R / (1 + 1i * k * R);

