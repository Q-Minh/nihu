function M = helmholtz_matrix_sp(I, J, c, k)
x = c(I,:);
y = c(J,:);
r = y-x;
r = sqrt(dot(r,r,2));
M = exp(1i*k*r)./r;
M(isinf(M)) = 0;
end
