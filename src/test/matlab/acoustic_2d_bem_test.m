clear
R0 = 1;
R1 = 2;
N = 100;

dir = [1 1 0];

radiator = create_circle_boundary(R0, N);
[x, n] = centnorm(radiator);
dir = dir / norm(dir);

field_mesh = create_circle_boundary(R1, N);

k = .5*mesh_kmax(radiator);

[r_nodes, r_elements] = extract_core_mesh(radiator);
[f_nodes, f_elements] = extract_core_mesh(field_mesh);

tic;
[Ls, Ms, Lf, Mf] = acoustic_2d_bem(r_nodes, r_elements, f_nodes, f_elements, k);
t_Booni = toc;

pinc = exp(-1i*k*(x*dir'));
qinc = -1i*k*(n*dir').*pinc;

qref = -qinc;
pref = Ms \ (Ls * qref);

ptot = pinc + pref;
