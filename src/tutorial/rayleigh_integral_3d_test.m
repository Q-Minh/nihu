clear;
clc;

surface = create_circle(1, 10);
field = create_slab([-2 0 .2; 2 0 .2; 2 0 2.2; -2 0 2.2], [40 20]);

[s_nodes, s_elem] = extract_Boonen_mesh(surface);
[f_nodes, f_elem] = extract_Boonen_mesh(field);

k = 5;

tic;
[Zf, Zs] = rayleigh_integral_3d(s_nodes, s_elem, f_nodes, f_elem, k);
toc;

vn = ones(size(Zf,2),1);    % // excitation
pf = Zf * vn;
ps = Zs * vn;

figure;
plot_mesh(surface, 20*log10(abs(ps)));
plot_mesh(field, 20*log10(abs(pf)));
