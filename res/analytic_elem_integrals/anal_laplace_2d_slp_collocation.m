clear all;

y = sym('y', 'real');
z = sym('z', 'positive');
d1 = sym('d1', 'positive');
d2 = sym('d2', 'positive');
eps = sym('eps', 'positive');

%% single layer potential
r = abs(y);
G = -log(r)/2/pi;

I = limit(int(G, y, -d1, -eps) + int(G, y, eps, d2), eps, 0);
pretty(I)


%% double layer potential
rvec = [y, -z];
r = sqrt(rvec * rvec.');
normal = [0 1];
G = -1/r^2 * rvec * normal.';

I = int(G, y, -d1, d2);


%% the hypersingular
r = abs(y);
G = 1/r^2/2/pi;
I = int(G, y, -d1, -eps) + int(G, y, eps, d2);
pretty(I)
