syms xi eta v

v = [1 xi xi^2 eta eta^2 xi*eta];

coords = [
    0 0
    .5 0
    1 0
    .5 .5
    0 1
    0 .5
    ];

for i = 1 : size(coords,1)
    V(i,:) = subs(subs(v, xi, coords(i,1)), eta, coords(i,2));
end

a = V \ eye(size(V));

N = simple(v * a).';

dN = simple([diff(N,xi), diff(N,eta)]);