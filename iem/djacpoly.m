function y = djacpoly(n,a,b,z,k)
% kth order derivative of Jacobian polynomial
if n == 0 
    y = zeros(size(z));
else
    y = gamma(a+b+n+1+k)/(2^k*gamma(a+b+n+1))*jacpoly(n-k,a+k,b+k,z);
end