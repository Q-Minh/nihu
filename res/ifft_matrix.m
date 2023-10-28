function M = ifft_matrix(L)
M = zeros(2*L+1, 2*(2*L)+1);
for i = 1 : 2*(2*L)+1
    x = zeros(2*(2*L)+1, 1);
    x(i) = 1;
    M(:,i) = padifft(x);
end
end
