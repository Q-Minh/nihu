function [p, roots, weights] = gauss_general(weight_fun, a, b, N)

syms x;

for n = 0 : N
    disp(n);
    p(n+1) = x^n;
    for k = 0 : n-1
        p(n+1) = p(n+1) - inner(p(k+1), p(n+1), weight_fun, a, b)*p(k+1);
    end
    p(n+1) = simple(normalise(p(n+1), weight_fun, a, b));
    
    if (n == 0)
        continue;
    end
    
    [r, L] = lagrange(p(n+1));
    w = int(L*weight_fun, a, b) ./ subs(weight_fun, x, r);

    [roots{n+1}, i] = sort(double(r));
    weights{n+1} = double(w(i));
end

end

function [r, L] = lagrange(p)
syms x;
r = solve(p);
total = prod(r-x);
L = total ./ (r-x);
L = L ./ diag(subs(L, x, r.'));
end

function I = inner(p, q, weight, a, b)
I = int(weight*p*q,'x',a,b);
end

function I = norm2(p, weight, a, b)
I = inner(p, p, weight, a, b);
end

function p = normalise(p, weight, a, b)
p = p / sqrt(norm2(p, weight, a, b));
end

