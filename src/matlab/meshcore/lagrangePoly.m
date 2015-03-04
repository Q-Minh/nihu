function [N, f, df] = lagrangePoly(base, corners)

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
N = s.' / q;

f = matlabFunction(N, 'vars', {g});
for j = 1 : d
    df{j} = matlabFunction(diff(N, g(j)), 'vars', {g});
end

end
