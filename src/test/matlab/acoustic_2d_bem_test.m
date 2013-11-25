clear

%% create meshes
R0 = 1;
R1 = 2;
N = 200;

radiator = create_circle_boundary(R0, N);
field_mesh = create_circle_boundary(R1, N);
[xs, ns] = centnorm(radiator);
[xf, nf] = centnorm(field_mesh);

[r_nodes, r_elements] = extract_core_mesh(radiator);
[f_nodes, f_elements] = extract_core_mesh(field_mesh);


%% wave number and system matrices
kmax = min(mesh_kmax(radiator));
k = .2 * kmax;

tic;
[Ls, Ms, Lf, Mf] = acoustic_2d_bem(r_nodes, r_elements, f_nodes, f_elements, k);
t_Booni = toc;

%% plane wave scattering

dir = [1 0 0];
dir = dir / norm(dir);

[pinc_s, qinc_s] = incident('plane', dir, xs, ns, k);
pinc_f = incident('plane', dir, xf, [], k);

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
xlabel('angle [rad]');
ylabel('scattered field');
legend('NiHu', 'Oltean');

%% 
x0 = [.2 -.1 0];

[ps_anal, qs] = incident('line', x0, xs, ns, k);
pf_anal = incident('line', x0, xf, [], k);

ps = Ms \ (Ls * qs);
pf = Mf * ps - Lf * qs;

figure;
plot(atan2(xs(:,2), xs(:,1)), abs(ps), '.', ...
    atan2(xs(:,2), xs(:,1)), abs(ps_anal), '.', ...
    atan2(xf(:,2), xf(:,1)), abs(pf), '.', ...
    atan2(xf(:,2), xf(:,1)), abs(pf_anal), '.');
xlabel('angle [rad]');
ylabel('scattered field');
legend('NiHu surface', 'Analytic', 'NiHu field', 'Analytic');

%% check matrices with matlab integration

[gx, nx, w, i] = geo2gauss(radiator, 3);
T = sparse(1:length(w), i, w);
[L, M] = incident('line', xf, gx, nx, k);
L = L.' * T;
M = M.' * T;

figure;

subplot(1,2,1);
pcolor(log10(abs(L./Lf-1)));
shading interp;
view(2);
colorbar;

subplot(1,2,2);
pcolor(log10(abs(M./Mf-1)));
shading interp;
view(2);
colorbar;
