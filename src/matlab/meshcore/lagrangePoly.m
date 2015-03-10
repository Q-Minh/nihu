function [N, dN, f, df, ddf] = lagrangePoly(base, corners)

d = size(corners,2);
n = size(corners,1);

g = sym('xi', [1 d]);
s = ones(n,1);
for j = 1 : d
    s = s .* g(j).^base(:,j);
end


q = ones(n,n);
for k = 1 : n
    for j = 1 : d
        q(:,k) = q(:,k) .* corners(:,j).^base(k,j);
    end
end
N = simplify(s.' / q);
for j = 1 : d
    dN(j,:) = simplify(diff(N, g(j)));
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
