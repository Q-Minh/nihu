clear;

surface = create_circle_boundary(.7, 100);
surface = flip_elements(surface);
field = create_circle_boundary(.3, 100);

x0 = [2 .3 0];

[c, n] = centnorm(surface);
rvec = bsxfun(@minus, c, x0);
r = sqrt(dot(rvec, rvec, 2));
ps_anal = -log(abs(r))/2/pi;
qs_anal = -1./r./r .* dot(rvec, n, 2);
[cf, nf] = centnorm(field);
rvec = bsxfun(@minus, cf, x0);
r = sqrt(dot(rvec, rvec, 2));
pf_anal = -log(abs(r))/2/pi;
qf_anal = -1./r./r .* dot(rvec, nf, 2);


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

