clear;
csound = 340 * (1+1e-2*1i);

%%
fs = 5000;
T = 10;
v = 30;    % source speed
t = -T/2:1/fs:T/2;  % time scale
t = t(1:end-1);
N = length(t);
y = 5;          % source y distance

f = (-N/2:N/2-1)/N*fs;
k = 2*pi*f/csound;

f0 = 1000;

%%
gtil = zeros(size(f));
kxx = 2*pi*(f-f0)/v;
arg = k.^2 - kxx.^2;
xpos = real(k).^2 - kxx.^2 > 0;
gtil(~xpos) = 1/2/pi * besselk(0, sqrt(-arg(~xpos))*y);
gtil(xpos) = -1i/4 * besselh(0, 2, sqrt(arg(xpos))*y);

xvec = [0 1 2 3]';
phat = repmat(gtil, length(xvec), 1) .* exp(-1i*xvec*kxx);

%%
p = fftshift(ifft(fftshift(phat,2), [], 2), 2);
