clear;
syms x y;
v = [1 x y x*x x*y y*y];

c = [
    0 0
    .5 0
    1 0
    .5 .5
    0 1
    0 .5
    ];

for i = 1 : size(c,1)
    V(i,:) = subs(subs(v, x, c(i,1)), y, c(i,2));
end
N = simple(v / V).';
dN = simple([diff(N,x) diff(N,y)]);
ddN = simple([diff(dN(:,1), x, 1), diff(dN(:,2), x, 1), diff(dN(:,2), y, 1)]);
