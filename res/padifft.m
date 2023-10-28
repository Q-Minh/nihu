function x = padifft(X)
    L = (length(X)-1) / 4;
    x_pad = ifft(X);
    x = zeros(2*L+1,1);
    x(1:L) = x_pad(end-L+1:end);
    x(L+1:end) = x_pad(1:L+1);
end
