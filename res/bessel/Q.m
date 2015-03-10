function p = Q(nu, z, N)
q = -1./(2*z).^2;
p = 0;
for m = 0 : N
    p = p + helper(nu,2*m+1) * q.^m;
end
p = p ./ (2*z);
end
