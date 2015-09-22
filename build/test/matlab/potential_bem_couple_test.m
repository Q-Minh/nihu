clear;

surface = create_sphere_boundary(1, 5);
field = create_sphere_boundary(3, 4);

[ns, es] = extract_core_mesh(surface);
[nf, ef] = extract_core_mesh(field);

tic;
[Gs, Hs, Gf, Hf] = potential_bem_couple(ns, es, nf, ef);
toc;

N = size(Gs,1);
qs = ones(N,1);
ps = Hs \ (Gs * qs);
pf = Hf * ps - Gf * qs;

clear mex
