clear;

surface = create_circle(1, 5);
surface = drop_unused_nodes(get_boundary(surface));
field = create_circle(3, 4);
field = drop_unused_nodes(get_boundary(field));

x0 = [.2 .3 .2];

[c, n] = centnorm(surface);
[ps_anal, qs_anal] = incident('point', x0, c, n, 0);
[cf, nf] = centnorm(field);
[pf_anal, qf_anal] = incident('point', x0, cf, nf, 0);

[nods, els] = extract_Boonen_mesh(surface);
[nodf, elf] = extract_Boonen_mesh(field);

tic;
[Gs, Hs, Gf, Hf] = potential_bem_2d(nods, els, nodf, elf);
toc;

ps = Hs \ (Gs * qs_anal);
err_ps = log10(abs(ps./ps_anal-1));

pf = Hf * ps - Gf * qs_anal;
err_pf = log10(abs(pf./pf_anal-1));

clear mex

