clear
R0 = 1;
R1 = 2;
N = 200;

dir = [1 1 0];
dir = dir / norm(dir);

radiator = create_circle_boundary(R0, N);
field_mesh = create_circle_boundary(R1, N);
[xs, ns] = centnorm(radiator);
[xf, nf] = centnorm(field_mesh);

k = 2;

[r_nodes, r_elements] = extract_core_mesh(radiator);
[f_nodes, f_elements] = extract_core_mesh(field_mesh);

tic;
[Ls, Ms, Lf, Mf] = acoustic_2d_bem(r_nodes, r_elements, f_nodes, f_elements, k);
t_Booni = toc;

pinc_s = exp(-1i*k*(xs*dir'));
qinc_s = -1i*k*(ns*dir').*pinc_s;
pinc_f = exp(-1i*k*(xf*dir'));

qref_s = -qinc_s;
pref_s = Ms \ (Ls * qref_s);
ptot_s = pinc_s + pref_s;

pref_f = Mf * pref_s - Lf * qref_s;
ptot_f = pinc_f + pref_f;

rf = sqrt(dot(xf,xf,2));
phif = atan2(xf(:,2), xf(:,1));
ptot_f_anal = planewave_cyl2d(rf, phif, R0, k, 30);

plot(atan2(xf(:,2), xf(:,1)), abs(ptot_s), '.', ...
    atan2(xf(:,2), xf(:,1)), abs(ptot_f_anal), '.');
