clear;
re = 0:1e-1:8;
im = -8:1e-1:1;
[Y,X] = meshgrid(im,re);
Z = X + 1i*Y;
tic;
[H0, H1] = bessel(Z);
tnum = toc;
tic;
H0anal = besselh(0,2,Z);
H1anal = besselh(1,2,Z);
tanal = toc;

figure;
subplot(2,3,1); pcolor(X,Y,real(H0)); shading flat;
subplot(2,3,2); pcolor(X,Y,imag(H0)); shading flat;
subplot(2,3,3); pcolor(X,Y,log10(abs(H0./H0anal-1))); shading flat; colorbar;
subplot(2,3,4); pcolor(X,Y,real(H1)); shading flat;
subplot(2,3,5); pcolor(X,Y,imag(H1)); shading flat;
subplot(2,3,6); pcolor(X,Y,log10(abs(H1./H1anal-1))); shading flat; colorbar;
