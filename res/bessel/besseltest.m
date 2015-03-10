clear all;
close all;
Z = 100;
x = linspace(-Z,Z,501);
nu = 0;
z = bsxfun(@plus, x, 1i*x');
F = BesselH(nu, 1, z);
f = besselh(nu, 1, z);

figure;
surf(x, x, log10(abs(F./f-1)));
view(2);
shading interp;
caxis([-15 -6]);
axis equal tight;
colorbar;
