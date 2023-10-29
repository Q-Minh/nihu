clear;
re = 0:1e-1:8;
im = -8:1e-1:-1e-3;
Z = bsxfun(@plus, re, 1i*im');
tic;
[J0, J1, Y0, Y1, H01, H11, H02, H12, K0, K1] = bessel(Z);
tnum = toc;
tic;
J0anal = besselj(0,Z);
J1anal = besselj(1,Z);
Y0anal = bessely(0,Z);
Y1anal = bessely(1,Z);
H01anal = besselh(0,1,Z);
H11anal = besselh(1,1,Z);
H02anal = besselh(0,2,Z);
H12anal = besselh(1,2,Z);
K0anal = besselk(0,Z);
K1anal = besselk(1,Z);
tanal = toc;

figure;
subplot(2,5,1); surf(re,im,log10(abs(J0./J0anal-1))); shading flat; caxis([-15 -6]); view(2); axis equal tight;
subplot(2,5,6); surf(re,im,log10(abs(J1./J1anal-1))); shading flat; caxis([-15 -6]); view(2); axis equal tight;
subplot(2,5,2); surf(re,im,log10(abs(Y0./Y0anal-1))); shading flat; caxis([-15 -6]); view(2); axis equal tight;
subplot(2,5,7); surf(re,im,log10(abs(Y1./Y1anal-1))); shading flat; caxis([-15 -6]); view(2); axis equal tight;
subplot(2,5,3); surf(re,im,log10(abs(H01./H01anal-1))); shading flat; caxis([-15 -6]); view(2); axis equal tight;
subplot(2,5,8); surf(re,im,log10(abs(H11./H11anal-1))); shading flat; caxis([-15 -6]); view(2); axis equal tight;
subplot(2,5,4); surf(re,im,log10(abs(H02./H02anal-1))); shading flat; caxis([-15 -6]); view(2); axis equal tight;
subplot(2,5,9); surf(re,im,log10(abs(H12./H12anal-1))); shading flat; caxis([-15 -6]); view(2); axis equal tight;
subplot(2,5,5); surf(re,im,log10(abs(K0./K0anal-1))); shading flat; caxis([-15 -6]); view(2); axis equal tight;
subplot(2,5,10); surf(re,im,log10(abs(K1./K1anal-1))); shading flat; caxis([-15 -6]); view(2); axis equal tight;
