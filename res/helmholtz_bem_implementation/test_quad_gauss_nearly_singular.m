clear;

%% Create element and field
coords = [
    0 0 0
    2 0 0
    2 2 0
    0 2 0
    ];
lset = ShapeSet.LinearQuad;

orders = [10 12];
zvec = .5;

integrals = nan(length(orders), length(zvec), 4);

for o = 1 : length(orders)
    order = orders(o);
    [xi, w] = gaussquad2(order, 4);
    [L, dL] = lset.eval(xi);
    x = L * coords;
    dxxi = dL(:,:,1) * coords;
    dxeta = dL(:,:,2) * coords;
    jvec = cross(dxxi, dxeta, 2);
    j = sqrt(dot(jvec, jvec, 2));
    ny = bsxfun(@times, jvec, 1./j);
    
    N = quad1_gauss_shapeset(xi);
    
    test_coords = coords;
    test_coords(:,3) = 1e-3;

    xi0vec = gaussquad2(2, 4);
    x0vec = lset.eval(xi0vec) * test_coords;
    N0 = quad1_gauss_shapeset(xi0vec);
    nx = [0 0 1];
    
    k = 1;
    K = @(x, y)laplace_3d_hsp_kernel(x, nx, y, ny);
    
    for i = 1 : length(zvec)
        progbar(1, length(zvec), i);
        x0 = x0vec(1,:);
        x0(3) = zvec(i);
        g = K(x0, x);
        
        integrand = bsxfun(@times, g .* j .* w, N);
        integrand0 = bsxfun(@times, g .* j .* w, N(1,:));
        integral(o,i,:) = sum(integrand - integrand0, 1);
    end
end
