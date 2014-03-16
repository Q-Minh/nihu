function K = elasto_3d_regular(type, mater, x0, coords)
switch type
    case 'T'
        kernel = @Tkernel;
    case 'U'
        kernel = @Ukernel;
end
n = cross(coords(2,:) - coords(1,:), coords(3,:) - coords(1,:));
jac = norm(n);
n = n / jac;

[xi, w] = gaussquad2(5, 3);
xg = shapefun(xi, 23) * coords;

K = zeros(3,3);
for g = 1 : size(xg,1)
    K = K + w(g) * kernel(mater, x0, xg(g,:), n);
end

K = jac * K;
end

