clear;
Le = 5e-2;
L = [1 1];
mesh = quad2tria(create_slab(L, ceil(L/Le)));
% mesh = quad2tria(create_circle(L(1), ceil(L(1)/Le)));
field_var = [2 1; 1 2];
space_var = [
    10 .1 0
    .1 .3 0
    0 0 1
    ];
[nodes, elements] = extract_core_mesh(mesh);
[D, B] = kltria_surface(nodes, elements, field_var, space_var);
D = (D + D') / 2;
B = (B + B') / 2;
nModes = 100;
[g, Lambda] = eigs(D, B, nModes);
lambda = diag(Lambda);

%%
figure;
for n = 1 : 8
    subplot(2,4,n);
    plot_mesh(mesh, g(2:2:end,n));
    shading interp;
    hold on;
    quiver(nodes(:,1), nodes(:,2), g(1:2:end,n), g(2:2:end,n), 'Color', 'black');
    axis equal tight off
end

%% generate random process
nReals = 8;
reals = g * diag(real(sqrt(lambda))) * randn(nModes, nReals);
figure;
for n = 1 : nReals
    subplot(2,4,n);
    p = reals(:,n);
    p = reshape(p, 2, [])';
    P = sqrt(dot(p,p,2));
    plot_mesh(mesh, P);
    shading interp;
    hold on;
    quiver(nodes(:,1), nodes(:,2), p(:,1), p(:,2), 'Color', 'black');
    axis equal tight off
    caxis([0 max(abs(P))]);
end
