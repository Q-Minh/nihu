close all
clear;

shapeset = ShapeSet.LinearQuad;
coords = [
    0 0 0
    1 0 0
    .8 1 0
    0 1 0
    ];

x = [.5 .5 .1];

order = 20;
[xi, w] = gaussquad2(order, 4);
[N, dN, ddN] = shapeset.eval(xi);

y = N * coords;
dyxi = dN(:,:,1) * coords;
dyeta = dN(:,:,2) * coords;
ddyxixi = ddN(:,:,1,1) * coords;
ddyetaeta = ddN(:,:,2,2) * coords;
ddyxieta = ddN(:,:,1,2) * coords;
Jvec = cross(dyxi, dyeta, 2);

xi0 = be_invmap(coords.', x, shapeset, 1e-8, 100);

[N0, dN0, ddN0] = shapeset.eval(xi0(1:2));
x0 = N0 * coords;
dyxi0 = dN0(:,:,1) * coords;
dyeta0 = dN0(:,:,2) * coords;
J0vec = cross(dyxi0, dyeta0, 2);
J0 = norm(J0vec);
d = (x-x0) / J0vec;

dxi = bsxfun(@minus, xi, xi0(1:2));
[theta, rho] = cart2pol(dxi(:,1), dxi(:,2));

figure;
plot3(y(:,1), y(:,2), y(:,3), 'b.', ...
    x(1), x(2), x(3), 'r.', ...
    x0(1), x0(2), x0(3), 'r*');
axis equal;

A1 = bsxfun(@times, dyxi, cos(theta)) + bsxfun(@times, dyeta, sin(theta));
A2 = bsxfun(@times, ddyxixi, cos(theta).^2)...
    + bsxfun(@times, ddyetaeta, sin(theta).^2) ...
    + 2 * bsxfun(@times, ddyxieta, sin(theta).*cos(theta));

rvec = bsxfun(@minus, y, x);
r = sqrt(dot(rvec, rvec, 2));
rappr = sqrt((d* J0)^2 + rho.^2 .* (dot(A1, A1, 2) - 2*d*(A2 * J0vec.')));

figure;
plot3(xi(:,1), xi(:,2), 1./r, '.',...
    xi(:,1), xi(:,2), 1./rappr, '.');
