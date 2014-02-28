function T = Tkernel(mater, x, y, n)
rvec = y(:) - x(:);
r = norm(rvec);
rdy = rvec / r;
rdn = dot(rdy, n(:), 1);
B = 1 / (8*pi*(1-mater.nu)*r^2);
T = -B * rdn * ((1-2*mater.nu)*eye(3) + 3*(rdy*rdy.'));
for i = 1 : 3
    for k = 1 : 3
        T(i,k) = T(i,k) + (1-2*mater.nu)*B * (rdy(k)*n(i) - rdy(i)*n(k));
    end
end
end
