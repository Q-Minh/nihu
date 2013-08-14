clear;

surface = create_sphere_boundary(1, 5);
surface = quad2tria(surface);
field = create_sphere_boundary(3, 4);

x0 = [0 0 0];

k = min(mesh_kmax(surface));

[cs, ns] = centnorm(surface);
[ps_anal, qs_anal] = incident('point', x0, cs, ns, k);
[cf, nf] = centnorm(field);
[pf_anal, qf_anal] = incident('point', x0, cf, nf, k);

[nods, els] = extract_Boonen_mesh(surface);
[nodf, elf] = extract_Boonen_mesh(field);

tic;
[Gs, Hs, Hts, Ks, Gf, Hf, Htf, Kf] = acoustic_bem(nods, els, nodf, elf, k);
toc;

ps = Hs \ (Gs * qs_anal);
err_ps = log10(abs(ps./ps_anal-1));

ps2 = Ks \ (Hts * qs_anal);
err_ps2 = log10(abs(ps2./ps_anal-1));

alph = 1i / k;
ps_bm = (Hs + alph * Ks) \ ((Gs + alph * Hts) * qs_anal);
err_ps_bm = log10(abs(ps_bm./ps_anal-1));

pf = Hf * ps - Gf * qs_anal;
err_pf = log10(abs(pf./pf_anal-1));

qf = Kf * ps - Htf * qs_anal;
err_qf = log10(abs(qf./qf_anal-1));

clear mex
