model = create_sphere_boundary(1, 10);
[nodes, element] = extract_Boonen_mesh(model);
tic;
M = mass_matrix(nodes, element);
toc;
fprintf('%g\n', sum(M(:))/4/pi);
clear mex
