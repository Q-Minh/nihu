mex -v CXXFLAGS="\$CXXFLAGS -std=c++0x" Boonen13.mex.cpp -I../../../eigen -output Boonen13

%%

mesh = create_sphere_boundary(1,5);
[nodes, elements] = extract_bem_mesh(mesh);
tic;
G = Boonen13(nodes, elements, 1);
tG = toc;
