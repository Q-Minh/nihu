clear;
syms x y;
v = x.^(0:2).' * y.^(0:2);
v = v(:);
v = v(1:end-1).';

c = [
    -1 -1
    0 -1
    1 -1
    1 0
    1 1
    0 1
    -1 1
    -1 0
    ];

for i = 1 : size(c,1)
    V(i,:) = subs(subs(v, x, c(i,1)), y, c(i,2));
end
N = simple(v / V).';
dN = simple([diff(N,x) diff(N,y)]);
