function T = Tkernel(mater, x, y, n)
rvec = y(:) - x(:);
r = norm(rvec);
rdy = rvec / r;
rdn = dot(rdy, n(:), 1);
T = -rdn * ((1-2*mater.nu)*eye(3) + 3*(rdy*rdy.')) +...
    (1-2*mater.nu)*(rdy*n.'-n*rdy.');
T = T / (8*pi*(1-mater.nu)*r^2);
end
