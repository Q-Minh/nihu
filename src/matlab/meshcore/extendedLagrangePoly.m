function [N, dN, f, df, ddf] = extendedLagrangePoly(base, corners)

N = lagrangePoly(base, corners);

for j = 1 : d
    dN(j,:) = simplify(diff(N, g(j)));
end

for j = 1 : d
    for k = 1 : d
        ddN(j,k,:) = simplify(diff(dN(j,:),g(k)));
    end
end

f = matlabFunction(N, 'vars', {g});
df = cell(n, d);
ddf = cell(n, d, d);
for i = 1 : n
    for j = 1 : d
        df{i,j} = matlabFunction(dN(j,i), 'vars', {g});
        for k = 1 : d
            ddf{i,j,k} = matlabFunction(simplify(diff(dN(j,i), g(k))), 'vars', {g});
        end
    end
end

end
