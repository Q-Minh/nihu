function [x, w] = laguerre(n)
%LAGUERRE Generate n-point Laguerre quadrature
%   LAGUERRE(N) prints the abscissae and weights of the N point Laguerre
%   quadrature that approximates int_0^1 f(x) dx where f(x) contains a
%   log(x) type singularity.

syms x;
L(0+1) = sym('1');
L(1+1) = 1-x;
for k = 1 : n
    L(k+1+1) = simple(1/(k+1)*((2*k+1-x)*L(k+1)-k*L(k-1+1)));
end

ti = solve(L(n+1));
w = ti/(n+1)/(n+1)./subs(L(n+1+1), x, ti).^2;
[x, i] = sort(double(exp(-ti)));
w = double(w(i));

if nargout == 0
    fprintf('\n');
    fprintf(1, '%.20x, %.20x,\n', [x w].');
end

end
