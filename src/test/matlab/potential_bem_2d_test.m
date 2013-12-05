clear;

surface = create_circle_boundary(1, 100);
field = create_circle_boundary(2, 100);

x0 = [.2 .3 0];

[c, n] = centnorm(surface);
rvec = bsxfun(@minus, c, x0);
r = sqrt(dot(rvec, rvec, 2));
ps_anal = -log(abs(r))/2/pi;
qs_anal = -1./r./r .* dot(rvec, n, 2) / 2/ pi;
[cf, nf] = centnorm(field);
rvec = bsxfun(@minus, cf, x0);
r = sqrt(dot(rvec, rvec, 2));
pf_anal = -log(abs(r))/2/pi;
qf_anal = -1./r./r .* dot(rvec, nf, 2) / 2/pi;

[nods, els] = extract_core_mesh(surface);
[nodf, elf] = extract_core_mesh(field);

tic;
[Ls, Ms, Mts, Ns, Lf, Mf] = potential_bem_2d(nods, els, nodf, elf);
toc;

ps = Ms \ (Ls * qs_anal);
err_ps = log10(abs(ps./ps_anal-1));

pf = Mf * ps - Lf * qs_anal;
err_pf = log10(abs(pf./pf_anal-1));

clear mex

