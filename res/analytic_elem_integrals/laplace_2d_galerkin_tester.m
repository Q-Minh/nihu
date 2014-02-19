function laplace_2d_galerkin_tester
x1 = [1 0];
x0 = [0 0];
x2 = [1 1];
[d1, d2, phi] = geom(x1, x0, x2);
fprintf(1, 'd1: %g, d2: %g, phi: %g\n', d1, d2, phi/pi*180);
I1 = anal_slp(d1, d2, phi);
I2 = anal_dlp(d1, d2, phi);
fprintf(1, 'I1: %g, I2: %g\n', I1, I2);
end

function [d1, d2, phi] = geom(x1, x0, x2)
v1 = x0-x1;
v2 = x2-x0;
d1 = norm(v1);
d2 = norm(v2);
phi = atan2(det([v1;v2]), v1 * v2.');
end

function I = anal_slp(d1, d2, phi)
c = cos(phi);
s = sin(phi);
d3 = sqrt(d1^2 + 2*d1*d2*c + d2^2);
I = (...
    d1*d2 * (3-2*log(d3)) ...
    + c*(d1^2*log(d1/d3) + d2^2*log(d2/d3))...
    - s*((d1^2*q(d2/d1,phi) + d2^2*q(d1/d2, phi)))...
    ) / (4*pi);
end

function I = anal_dlp(d1, d2, phi)
c = cos(phi);
s = sin(phi);
d3 = sqrt(d1^2 + 2*d1*d2*c + d2^2);
I = (...
    d2*c * q(d1/d2, phi)...
    -d1 * q(d2/d1, phi)...
    - d2*s*log(d3./d2))/(2*pi);
end

function f = q(a, phi)
f = atan(a/sin(phi)+cot(phi))-atan(cot(phi));
end
