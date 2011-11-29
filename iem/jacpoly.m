function y = jacpoly(n,a,b,z)
%JACPOLY Calculates Jacobi polynomials of order n, type (a,b) at
%coordinates z. z must be a column vector.
m = 0:n;
nz = length(z);
temp = repmat(binomial(n,m).*factorial(a+b+n+m)./factorial(a+m),nz,1).*repmat((z-1)/2,1,n+1).^repmat(m,nz,1);
y = repmat(gamma(a+n+1)/(factorial(n)*gamma(a+b+n+1)),nz,1).*sum(temp,2);