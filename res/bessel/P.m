function p = P(nu, z, N)
q = -1/(2*z).^2;
p = 0;
for m = 0 : N
    p = p + helper(nu,2*m) * q.^m;
end
end
