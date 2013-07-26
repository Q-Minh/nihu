clear;

surface = create_sphere_boundary(1, 5);
field = create_sphere_boundary(3, 4);

[ns, es] = extract_Boonen_mesh(surface);
[nf, ef] = extract_Boonen_mesh(field);

tic;
[Gs, Hs, Hs2, Gf, Hf, Hf2] = potential_bem_triple(ns, es, nf, ef);
toc;

N = size(Gs,1);
qs = ones(N,1);
ps = Hs \ (Gs * qs);
pf = Hf * ps - Gf * qs;

clear mex

