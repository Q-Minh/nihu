clear;
%close all;

tria = ~true;
gmsh = true;

Lout = [.2 .2 1];
Lin = [.1 .1 .9];
Le = 1e-1/2;
eps = 1e-3;

if ~gmsh
mesh = create_brick(Lout, ceil(Lout/Le));
mesh = translate_mesh(mesh, -Lout/2);

mesh = drop_mesh_IDs(drop_unused_nodes(get_boundary(mesh)));
mesh = flip_mesh(mesh);

if (tria)
    mesh = quad2tria(mesh);
end
else
    mesh = import_gmsh_mesh('pipe_03.msh');
    sel = mesh.Elements(:,2) == 23;
    mesh.Elements(sel, 2) = ShapeSet.LinearTria.Id;
    sel = mesh.Elements(:,2) == 24;
    mesh.Elements(sel, 2) = ShapeSet.LinearQuad.Id;
    mesh = flip_mesh(mesh);
end

c = 340;
kmax = min(mesh_kmax(mesh, 8));

fmax = kmax * c / (2*pi);
fprintf('Max mesh freq: %.1f Hz\n', fmax);

[r_nodes, r_elems] = extract_core_mesh(mesh);

field_xy = translate_mesh(create_slab(Lin([1 2]), Lin([1 2])/Le), -Lin([1 2])/2);
field_xz = translate_mesh(create_slab(Lin([1 3]), Lin([1 3])/Le), -Lin([1 3])/2);
field_xz = rotate_mesh(field_xz, pi/2, [1 0 0]);
field_yz = translate_mesh(create_slab(Lin([3 2]), Lin([3 2])/Le), -Lin([3 2])/2);
field_yz = rotate_mesh(field_yz, pi/2, [0 1 0]);

field = join_meshes(field_xy, field_xz, field_yz);
[f_nodes, f_elems] = extract_core_mesh(field);

[~, exc_sel] = mesh_select(mesh, sprintf('abs(z) > %f', Lout(3)/2-eps), 'ind', 'all');
%%


%// call C++ code at wave number k \approx pi
fvec = 10 : 5 : 600;
nFreq = length(fvec);


nElem = size(mesh.Elements, 1);
nField = size(field.Elements, 1);

ps_conv = zeros(nElem, nFreq);
pf_conv = zeros(nField, nFreq);

ps_bm = zeros(nElem, nFreq);
pf_bm = zeros(nField, nFreq);

ps_hs = zeros(nElem, nFreq);
pf_hs = zeros(nField, nFreq);

qs_ana = zeros(nElem, 1);
qs_ana(exc_sel) = 1;

ps_ana = zeros(nElem, nFreq);
pf_ana = zeros(nField, nFreq);

r_c = centnorm(mesh);
f_c = centnorm(field);

err_s_conv = zeros(1, nFreq);
err_f_conv = zeros(1, nFreq);

err_s_hs = zeros(1, nFreq);
err_f_hs = zeros(1, nFreq);

err_s_bm = zeros(1, nFreq);
err_f_bm = zeros(1, nFreq);


for iFreq = 1 : nFreq

    progbar(1, nFreq, iFreq);
    f = fvec(iFreq);
    k = 2*pi*f/c;
    [Ls, Ms, Mts, Ns, Lf, Mf] =...
        singular_bem_3d_field(r_nodes, r_elems, f_nodes, f_elems, k);

    %// acoustic pressure on the surface and in the field points
    ps_conv(:, iFreq) = Ms \ (Ls * qs_ana);               %// conventional
    coup = 1i/k;   %// coupling constant
    ps_bm(:, iFreq) = (Ms + coup * Ns) \ (Ls * qs_ana + coup * Mts * qs_ana);

    ps_hs(:, iFreq) = (Ns) \ (Mts * qs_ana);
    
    % Compute field
    pf_conv(:, iFreq) = Mf * ps_conv(:, iFreq) - Lf * qs_ana;
    pf_bm(:, iFreq) = Mf * ps_bm(:, iFreq) - Lf * qs_ana;
    pf_hs(:, iFreq) = Mf * ps_hs(:, iFreq) - Lf * qs_ana;
    
    % Solve analytically
    A = [ 1i*k*exp(-1i*k*Lout(3)/2), -(1i*k)*exp(1i*k*Lout(3)/2);
         -1i*k*exp(-1i*k*(-Lout(3)/2)), (1i*k)*exp(1i*k*(-Lout(3)/2)) ];
    sol = A \ [1; 1]; 
    pp = sol(1); pm = sol(2);
    
    ps_ana(:, iFreq) = pp*exp(-1i*k*r_c(:,3)) + pm*exp(1i*k*r_c(:,3));
    pf_ana(:, iFreq) = pp*exp(-1i*k*f_c(:,3)) + pm*exp(1i*k*f_c(:,3));
    
    err_s_conv(:, iFreq) = norm(ps_conv(:, iFreq) - ps_ana(:, iFreq)) / norm(ps_ana(:, iFreq));
    err_s_hs(:, iFreq) = norm(ps_hs(:, iFreq) - ps_ana(:, iFreq)) / norm(ps_ana(:, iFreq));
    err_s_bm(:, iFreq) = norm(ps_bm(:, iFreq) - ps_ana(:, iFreq)) / norm(ps_ana(:, iFreq));
    
    err_f_conv(:, iFreq) = norm(pf_conv(:, iFreq) - pf_ana(:, iFreq)) / norm(pf_ana(:, iFreq));
    err_f_hs(:, iFreq) = norm(pf_hs(:, iFreq) - pf_ana(:, iFreq)) / norm(pf_ana(:, iFreq));
    err_f_bm(:, iFreq) = norm(pf_bm(:, iFreq) - pf_ana(:, iFreq)) / norm(pf_ana(:, iFreq));
end
%%
if (0)
 figure;
 subplot(3, 1, 1);
  plot(fvec, real(pf_conv(1, :)), ...
     fvec, real(pf_bm(1, :)), ...
     fvec, real(pf_hs(1, :)));
 hold on;
 plot(fvec, real(pf_ana(1,:)), 'k');
 subplot(3, 1, 2);
 plot(fvec, imag(pf_conv(1, :)), ...
     fvec, imag(pf_bm(1, :)), ...
     fvec, imag(pf_hs(1, :)));
 hold on;
 plot(fvec, imag(pf_ana(1,:)), 'k');
subplot(3, 1, 3);
 plot(fvec, log10(abs(pf_conv(1, :))), ...
     fvec, log10(abs(pf_bm(1, :))), ...
     fvec, log10(abs(pf_hs(1, :))));
 hold on;
 plot(fvec, log10(abs(pf_ana(1,:))), 'k');
% // plot results
end
%%
figure;
semilogy(fvec, err_s_conv, fvec, err_s_bm, fvec, err_s_hs, 'LineWidth', 1.2);
legend({'Conventional', 'Burton-Miller', 'Hypersingular'});
xlabel('Frequency [Hz]');
ylabel('Relative error');
%%
if (1)
f_sel = 200;
f_idx = find(fvec == f_sel, 1, 'first');
figure;
subplot(2,2,1);
plot_mesh(mesh, real(ps_conv(:,f_idx)));
axis equal;
colorbar;
caxis([-1 1]*max(abs(ps_ana(:, f_idx))));
title('Conventional');
subplot(2,2,2);
plot_mesh(mesh, real(ps_bm(:,f_idx)));
axis equal;
colorbar;
caxis([-1 1]*max(abs(ps_ana(:, f_idx))));
title('Burton');
subplot(2,2,3);
plot_mesh(mesh, real(ps_hs(:,f_idx)));
axis equal;
colorbar;
caxis([-1 1]*max(abs(ps_ana(:, f_idx))));
title('Hyper');
subplot(2,2,4);
plot_mesh(mesh, real(ps_ana(:,f_idx)));
axis equal;
colorbar;
caxis([-1 1]*max(abs(ps_ana(:, f_idx))));
title(sprintf('Analytical f = %g Hz', f_sel));
end

%%
figure
plot(r_c(:, 3), real(ps_hs(:, f_idx)), '*'); hold on;
plot(r_c(:, 3), real(ps_ana(:, f_idx)), 'r*');
grid on;