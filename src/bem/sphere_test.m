%% compile C++ mex code
clc;
mex -v CXXFLAGS="\$CXXFLAGS -std=c++11 -O3" Boonen13.mex.cpp -I../../../eigen -output Boonen13

%% Build a good little mesh
R = 1;                                  % radius
mesh = create_sphere_boundary(1, 4);    % mesh
k = min(mesh_kmax(mesh)) / 2;           % wave number

%% C++ Boonen version
% qmesh = quadratise(mesh);
% qmesh.Nodes(:,2:4) = bsxfun(@times, qmesh.Nodes(:,2:4), 1./sqrt(dot(qmesh.Nodes(:,2:4), qmesh.Nodes(:,2:4), 2)));
tic;
[nodes, elements] = extract_Boonen_mesh(mesh);
[G, H] = Boonen13(nodes, elements, k);
tBooni = toc;

%% old C++ version
tic;
[nodes, elements] = extract_bem_mesh(mesh);
[HH, GG] = bemHG(mesh, k, 'const');
tOldSchool = toc;

%% Solve
M = ibemRHS(mesh);
p = (H - .5*M) \ sum(G,2);
pp = (HH - .5*eye(size(GG,1))) \ sum(GG,2);

%% Analytical
p_anal = -R/(1 + 1i*k*R);

