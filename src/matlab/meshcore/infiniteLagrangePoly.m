function [N, g] = infiniteLagrangePoly(base, corners)

d = size(corners,2);
n = size(corners,1);

g = sym('xi', [1 d]);

if (d == 1)
    N = sym(1);
else
    N = lagrangePoly(base, corners(1:n/2,1:d-1));
end

xi = g(d);
if d == 2
    N = simplify([N * (2*xi/(xi-1)) fliplr(N) * (1+xi)/(1-xi)]);
else
    N = simplify([N * (2*xi/(xi-1)) N * (1+xi)/(1-xi)]);
end

end
