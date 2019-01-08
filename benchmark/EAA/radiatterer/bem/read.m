clear;

%%
c = 340;
fmax = min(mesh_kmax(mesh)) * c / 2/pi;

%%
res = importdata('data/response.res');
k = res(1);
ps = complex(res(2:2:end), res(3:2:end));
ps = reshape(ps, N, []).';
figure;
plot_mesh(mesh, mean(abs(ps), 2));

