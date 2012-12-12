clear;

%%
ns = 1000;
nr = 500;
Dy = 4;
xs = rand(ns,3)-.5; % source location
Xs = [0 0 0];       % source multipole location
xr = rand(nr,3)-.5 + repmat([0 Dy 0], nr, 1);       % receiver multipole location
Xr = [0 Dy 0];    % receiver location
qs = rand(ns, 1);    % source strengths

f0 = 50;
csound = 340;
v = 30;

%%
fs = 1000;
T = 5;
t = -T : 1/fs : T;
t = t(1:end-1);
x = t*v;
N = length(t);
fvec = (-N/2 : N/2-1)/N * fs;
kx = 2*pi*fvec/v;

%%
ds = repmat(Xs, ns, 1) - xs;       % source local distance vector
dr = xr - repmat(Xr, nr, 1);       % receiver local distance vector
D = Xr - Xs;        % multipole distance vector

%%
Dclus = 1;
C = 3;
Dvec(:,1) = x;
Dvec(:,2) = Dy;
Dvec(:,3) = 0;

klim = abs(min(kx)*2);

g_multi = zeros(nr, N);
%%
figure;
h_plot = plot(fvec, 20*log10(abs(g_multi(1,:))));
%%
for iF = length(fvec)/2 + [100:500]
    f = fvec(iF);
    
    k = 2*pi*f/csound;
    L = round(abs(k*Dclus) + C*log10(abs(k*Dclus)+pi))+2;

    [S, w, ~, sphere] = spherequad(L);  % quadrature on unit sphere

    kxx = 2*pi*(f-f0)/v;
    
    Fs = (qs .* exp(1i*kxx*ds(:,1))).' * exp(-1i*k*(ds*S));         % far field signature
    M = translation(k, Dvec, L, S);    % translation
    Mtil = fftshift(ifft(fftshift(M,2), [], 2), 2);
    mtil = interp1(kx, Mtil.', mod(kxx+klim/2, klim)-klim/2);
    Ns = Fs .* mtil;
    gr = exp(-1i*k*(dr*S));         % receiver term
    g_multi(:,iF) = (gr .* (exp(1i*kxx*dr(:,1)) * Ns)) * w;    % total transfer
    
    set(h_plot, 'yData', 20*log10(abs(g_multi(1,:))));
    drawnow;
end

%%

p = fftshift(ifft(fftshift(g_multi, 2), [], 2), 2);

%%
g_anal = zeros(nr,1);
for i = 1 : ns
    g_anal = g_anal + qs(i) * incident('point', xs(i,:), xr, [], k);
end

error = abs(g_multi-g_anal)./abs(g_anal);
fprintf(1, 'Error: %x\n', error);