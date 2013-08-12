clear;

surface = create_sphere_boundary(1, 5);
surface = quad2tria(surface);
field = create_sphere_boundary(3, 4);

x0 = [.2 .3 .2];

[c, n] = centnorm(surface);
[ps_anal, qs_anal] = incident('point', x0, c, n, 0);
[cf, nf] = centnorm(field);
[pf_anal, qf_anal] = incident('point', x0, cf, nf, 0);

[nods, els] = extract_Boonen_mesh(surface);
[nodf, elf] = extract_Boonen_mesh(field);

tic;
[Gs, Hs, Hts, Ks, Gf, Hf, Htf, Kf] = potential_bem(nods, els, nodf, elf);
toc;

ps = Hs \ (Gs * qs_anal);
err_ps = log10(abs(ps./ps_anal-1));

ps2 = Ks \ (Hts * qs_anal);
err_ps2 = log10(abs(ps2./ps_anal-1));

pf = Hf * ps - Gf * qs_anal;
err_pf = log10(abs(pf./pf_anal-1));

qf = Kf * ps - Htf * qs_anal;
err_qf = log10(abs(qf./qf_anal-1));

clear mex
