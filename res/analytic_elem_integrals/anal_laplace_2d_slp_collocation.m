clear all;

y = sym('y');
R1 = sym('R1', 'positive');
R2 = sym('R2', 'positive');
eps = sym('eps', 'positive');

%% analytical derivation
r = abs(y);
G = -log(r)/2/pi;

I = limit(int(G, y, -R1, -eps) + int(G, y, eps, R2), eps, 0);
pretty(I)
