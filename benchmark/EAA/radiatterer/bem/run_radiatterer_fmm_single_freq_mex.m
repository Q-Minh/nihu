clear;

type = 'gauss';

le = 5e-2;
radiator = create_radiatterer(le);
f = 346.5;
c = 340;
rho = 1.3;
om = 2*pi*f;
k = om/c;

switch type
    case 'const'
        mex_fun = @helmholtz_3d_hffmm_mex;
        order = 1;
    case 'gauss'
        mex_fun = @helmholtz_3d_hffmm_gauss_mex;
        order = 2;
end

%% Initialize, setup parameters
mex_fun('init');

accuracy = 3.0;
mex_fun('set', 'accuracy', 3.0, 'wave_number', k);

%% Setup mesh
[r_nodes, r_elems] = extract_core_mesh(radiator, 'surface');
mex_fun('mesh', r_nodes, r_elems);


%% Setup tree
leaf_diameter = 1 / k;
mex_fun('tree', 'divide_num_nodes', 10);
%% Assemble matrices
fprintf('Assembling FMM matrices ... '); tic;
mex_fun('matrix');
fprintf('Ready in %.3f seconds\n', toc);

%% Prepare excitation and right hand side vector
[gc, gn] = geo2gauss(radiator, order);
q = -1i*om*rho*1e-3 * ones(size(gc,1),1);
rhs = mex_fun('mvp_slp', q);

%% Prepare solver 
fprintf('Solving FMBEM system ... '); tic;
Afun = @(x)mex_fun('mvp_dlp', ensure_complex(x));
[p_bem, flag, relres, iter, resvec] = gmres(Afun, rhs, [], 1e-8, 1000);
fprintf('Ready in %.3f seconds.\n', toc);

%% Evaluate and plot results
p_bem_plot = mean(reshape(p_bem, order^2, []), 1).';
figure;
plot_mesh(radiator, 20*log10(abs(p_bem_plot)/2e-5));
view(150,25);
shading flat;
axis equal;
hcb = colorbar;
ylabel(hcb, 'Sound pressure level [dBSPL]');
caxis([-40 0] + 110);
hl = light;
lighting phong;
hl.Position = [2 2 1];
%%
figure;
plot_mesh(radiator, abs(p_bem_plot));

%%
mex_fun('cleanup');

%%
%clear mex;
save(sprintf('data/%s_%03dmm/%s_%03dmm_%dHz_ps', ...
    type, round(le*1000), type, round(le*1000), round(f)), 'p_bem');