%// create meshes
R0 = 1;			%// radius of the cylinder
R1 = 2;			%// radius of the field point mesh
N = 200;		%// number of boundary elements on the cylinder
dir = [1 0 0];	%// direction of the incident plane wave

radiator = create_circle_boundary(R0, N);
field_mesh = create_circle_boundary(R1, N);

kmax = min(mesh_kmax(radiator));
k = .5 * kmax;

%// system matrices
[r_nodes, r_elements] = extract_core_mesh(radiator);
[f_nodes, f_elements] = extract_core_mesh(field_mesh);
[Ls, Ms, Lf, Mf] = helmholtz_bem_2d_cyl(...
    r_nodes, r_elements, f_nodes, f_elements, k);

%// Numerical solution
[xs, ns] = centnorm(radiator);	%// center and normal of the cylinder elements
[pinc_s, qinc_s] = incident('plane', dir, xs, ns, k);	%// incident wave field

qref_s = -qinc_s;				%// reflected velocity field
pref_s = Ms \ (Ls * qref_s);	%// solution of the BEM problem
ptot_s = pinc_s + pref_s;		%// total pressure field

[xf, nf] = centnorm(field_mesh);%// center and normal of the field elements
pinc_f = incident('plane', dir, xf, [], k);	%// incident wave field

pref_f = Mf * pref_s - Lf * qref_s;	%// radiated pressure to field points
ptot_f = pinc_f + pref_f;			%// total pressure in field points

%// Analytical solution
pref_f_anal = planewave_cyl2d(xf, R0, k, 100);
ptot_f_anal = pinc_f + pref_f_anal;

plot(atan2(xf(:,2), xf(:,1)), abs(ptot_f), '.-', ...
    atan2(xf(:,2), xf(:,1)), abs(ptot_f_anal), '.-');
xlabel('angle [rad]');
ylabel('scattered field');
legend('NiHu', 'Analytic');

error = log10(mean(abs(ptot_f./ptot_f_anal-1)));
fprintf(1, 'Log10 mean error: %f\n', error);
