function X = padfft(x)
    L = (length(x)-1) / 2;
    x_pad = zeros(2 * (2*L) + 1, 1);
    x_pad(1:L+1) = x(end-L : end);
    x_pad(end-L+1:end) = x(1:L);
    X = fft(x_pad);
end
