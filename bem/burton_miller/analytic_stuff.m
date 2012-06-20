syms k r;

Gk = exp(-1i*k*r) / (4*pi*r);

dGk = diff(Gk, 'r', 1);
ddGk = diff(Gk, 'r', 2);