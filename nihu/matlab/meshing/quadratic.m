function x = quadratic(a, b, c)
%QUADRATIC  Solve quadratic polynomial equation ax^2 + bx + c = 0

d = max(abs([a,b,c]));
if abs(a)/d > 1e-8
    det = sqrt(b^2-4*a*c);
    x(1) = (-b + det) / (2*a);
    x(2) = (-b - det) / (2*a);
else
    x(1) = -c/b;
    x(2) = Inf;
end
