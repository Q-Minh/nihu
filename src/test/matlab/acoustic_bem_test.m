clear;

surface = create_sphere_boundary(1, 5);
surface = quad2tria(surface);
field = create_sphere_boundary(3, 4);

[ns, es] = extract_Boonen_mesh(surface);
[nf, ef] = extract_Boonen_mesh(field);

tic;
[Gs, Hs, Gf, Hf] = acoustic_bem(ns, es, nf, ef);
toc;

N = size(Gs,1);
qs = ones(N,1);
ps = Hs \ (Gs * qs);
pf = Hf * ps - Gf * qs;

figure;
plot_mesh(surface, abs(ps));
plot_mesh(field, abs(pf));
alpha .5;

clear mex
