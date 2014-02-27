clear;

x = sym('x', [3, 1]);
y = sym('y', [3, 1]);
rvec = y - x;
r = sqrt(rvec.' * rvec);
rdy = rvec / r;
n = sym('n', [3 1]);

nu = sym('nu', 'positive');

% mu = sym('mu', 'positive');
% A = 1/16*pi*mu*(1-nu);
% B = (1-2*nu) / (8*pi*mu*(1-nu));
A = sym('A', 'positive');
B = sym('B', 'positive');

U = A/r * ((3-4*nu)*eye(3) + rdy * rdy.');
for i = 1 : 3
    for k = 1 : 3
        T(i,k) = B/r^2 * (rdy(k)*n(i) - rdy(i)*n(k));
    end
end
% U = 1/(A*r) * ((3-4*mu)*eye(3));

U = simple(subs(U, x, zeros(3,1)));
T = simple(subs(T, x, zeros(3,1)));

rho = sym('rho', 'positive');
theta = sym('theta', 'real');

U = simple(subs(U, y, [rho*cos(theta); rho*sin(theta); 0]));
T = simple(subs(T, y, [rho*cos(theta); rho*sin(theta); 0]));
T = simple(subs(T, n, [0; 0; 1]));

R = sym('R', 'positive');
eps = sym('eps', 'positive');

IU1 = simple(int(U * rho, rho, 0, R/cos(theta)));
IU2 = int(IU1, theta);

IT1 = simple(int(T * rho));
IT1 = subs(IT1, rho, R/cos(theta)) - subs(IT1, rho, eps);
IT1 = limit(IT1, eps, Inf); % cancel epsilon member
IT2 = int(IT1, theta);

pretty(simple(IU2/A/R))
pretty(simple(IT2/B*R))

