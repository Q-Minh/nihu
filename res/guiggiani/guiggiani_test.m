clear;

%%
X = [
    -1 0 1
    +1 0 0
    +1 +1 0
    -1 +1 0
    ];
xi0 = [0 0];

id = 24;
theta = 0.0;

[Fm2, Fm1] = guigginai_series_expansion(theta, xi0, id, X);
fprintf(1, '4 pi * F-2: %g\n', Fm2*4*pi);
fprintf(1, '4 pi * F-1: %g\n', Fm1*4*pi);

%%
Ianal = -2*sin(pi/4)*4/4/pi;

nvec = 2 : 10;
for i = 1 : length(nvec)
    I(i,:) = Guiggiani(nvec(i), 24, X, xi0);
end

plot(nvec, log10(abs(I/Ianal - 1)), '.-');
xlabel('n');
ylabel('log10 error');
