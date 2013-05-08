syms xi eta v

v = [1 xi xi^2 eta eta^2 xi*eta xi*eta^2 eta*xi^2 eta^2*xi^2];

coords = [
    -1 -1
    0 -1
    1 -1
    1 0
    1 1
    0 1
    -1 1
    -1 0
    0 0
    ];

for i = 1 : size(coords,1)
    V(i,:) = subs(subs(v, xi, coords(i,1)), eta, coords(i,2));
end

a = V \ eye(size(V));

N = simple(v * a).';

dN = simple([diff(N,xi), diff(N,eta)]);