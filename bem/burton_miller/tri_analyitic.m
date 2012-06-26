syms p t;

% Term J1
I1 = int(int(sin(t)*sin(t)*cos(t)*cos(p), 't',0, pi/2), 'p', 0, 2*pi);
I2 = int(int(sin(t)*sin(t)*cos(t)*sin(p), 't', 0, pi/2), 'p', 0, 2*pi);
I3 = int(int(sin(t)*cos(t)*cos(t), 't', 0, pi/2), 'p', 0, 2*pi);

syms eps k;

ee = exp(-1i*k*eps)/(2*eps);

l1 = limit(diff(ee,'eps',1)*eps, 'eps', 0);

% Term J22
I4 = int(int(sin(t)*cos(t), 't', 0, pi/2), 'p', 0, 2*pi);
