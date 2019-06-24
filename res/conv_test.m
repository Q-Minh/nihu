clear;
L = 15;
s = rand(2*L+1,1) + 1i * rand(2*L+1,1);
m = rand(2*L+1,1) + 1i * rand(2*L+1,1);

%%
l = conv(m,s);
l = l(2*L+1 + (-L : L));

S = fft(s, 2*(2*L)+1);
M = fft(m, 2*(2*L) + 1);
l_fft = ifft(S .* M);
l_fft = l_fft(2*L+1 + (-L : L));

S = padfft(s);
M = padfft(m);
l_pfft = padifft(S .* M);

plot(real([l l_fft l_pfft]));

%%
l = conv(m,flipud(s));
l = l(2*L+1 + (-L : L));

S = padfft(s);
M = conj(padfft(conj(m)));
l_pfft = padifft(S .* M);

figure;
subplot(2,1,1);
plot(real([l l_pfft]));
subplot(2,1,2);
plot(imag([l l_pfft]));

