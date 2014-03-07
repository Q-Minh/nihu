function U = Ukernel(mater, x, y, ~)
rvec = y(:) - x(:);
r = norm(rvec);
rdy = rvec / r;
A = 1/(16*pi*mater.mu*(1-mater.nu)*r);
U = A * ((3-4*mater.nu)*eye(3) + (rdy*rdy.'));
end
