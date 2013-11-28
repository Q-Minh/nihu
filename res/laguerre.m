clear;

n = 10;


syms x;
L(0+1) = sym('1');
L(1+1) = 1-x;
for k = 1 : n
    L(k+1+1) = simple(1/(k+1)*((2*k+1-x)*L(k+1)-k*L(k-1+1)));
end

ti = solve(L(n+1));
wi = ti/(n+1)/(n+1)./subs(L(n+1+1), x, ti).^2;
[xi, i] = sort(double(exp(-ti)));
wi = double(wi(i));

fprintf('\n');
fprintf(1, '%.20x, %.20x,\n', [xi wi].')