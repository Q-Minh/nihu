clear;

%%
X = [
    -1 -1 0
    +1 -1 0
    +1 +1 0
    -1 +1 0
    ];
xi0 = [0 0];

theta = 0.0;
shape = ShapeSet.LinearQuad;
[Fm2, Fm1] = guiggiani_series_expansion_3d(theta, xi0, shape, X);
fprintf(1, '4 pi * F-2: %g\n', Fm2*4*pi);
fprintf(1, '4 pi * F-1: %g\n', Fm1*4*pi);

%%
Ianal = -2*sin(pi/4)*4/4/pi;

nvec = 2:10;
for i = 1 : length(nvec)
    [I(i,:), Isurf, Ilin1, Ilin2] = Guiggiani_3d(nvec(i), shape, X, xi0);
    fprintf(1, 'I0: %g\nI1: %g\nI2: %g\nI: %g\n', Isurf, Ilin1, Ilin2, I);
end

plot(nvec, log10(abs(sum(I,2)/Ianal - 1)), '.-');
xlabel('n');
ylabel('log10 error');
