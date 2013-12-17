clear;

X = [
    -1 -1 0
    1 -1 0
    1 1 0
    -1 1 0
    ];
xi0 = [0 0];

Ianal = -2*sin(pi/4)*4/4/pi;

nvec = 2 : 10;
for i = 1 : length(nvec)
    I(i,:) = Guiggiani(nvec(i), 24, X, xi0);
end

plot(nvec, log10(abs(I/ Ianal - 1)), '.-');
xlabel('n');
ylabel('log10 error');
