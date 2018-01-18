clear;

xi = sym('xi', 'real');
eta = sym('eta', 'real');
nx = sym('nx', [1 2], 'real');
d = sym('d', 'positive');

rvec = [xi, -eta];
r = sqrt(dot(rvec, rvec, 2));
ny = [0 1];
rdny = dot(rvec, ny, 2) / r;
rdnx = -dot(rvec, nx, 2) / r;

G = -log(r);
Gy = -1/r * rdny;
Gx = -1/r * rdnx;
Gxy = 1/r^2*(dot(nx, ny, 2) + 2 * rdnx * rdny);
Gxx = 1/r^2 * (-1 + 2*rdnx^2);
Gxxy = 2/r^3 * (rdny - 2*rdnx*dot(nx, ny, 2) - 4*rdnx^2*rdny);

prim = @(f)simplify(int(f, xi));

IG = prim(G);
IGy = prim(Gy);
IGx = prim(Gx);
IGxy = prim(Gxy);
IGxx = prim(Gxx);
IGxxy = prim(Gxxy);

sing = @(f)limit(subs(subs(f, xi, d) - subs(f, xi, -d), nx, [0 1]), eta, 0);

IG0 = sing(IG);
IGy0 = sing(IGy);
IGx0 = sing(IGy);
IGxy0 = sing(IGxy);
IGxx0 = sing(IGxx);
IGxxy0 = sing(IGxxy);