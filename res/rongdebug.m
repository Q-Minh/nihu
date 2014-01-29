clear;

order = 9;
nGauss = floor(order/2+1);

xi = sym('xi', [2 1]);
eta = sym('eta', [2 1]);

Corners = [
    -1 1 1 -1
    -1 -1 1 1
    ];

L = [
    (1-xi(1))*(1-xi(2)), ...
    (1+xi(1))*(1-xi(2)), ...
    (1+xi(1))*(1+xi(2)), ...
    (1-xi(1))*(1+xi(2))
    ] / 4;
dL = [
    diff(L, xi(1), 1)
    diff(L, xi(2), 1)
    ];
ddL = [
    diff(dL(1,:), xi(1), 1)
    diff(dL(2,:), xi(1), 1)
    diff(dL(2,:), xi(2), 1)
    ];

X = [
    0 .5 1 0
    0 0 1 1
    0 0 0 0
    ];

xi0 = [0 0].';

x0 = X * subs(L, xi, xi0).';
dx0 = X * subs(dL, xi, xi0).';

gamma = acos(dot(dx0(:,1), dx0(:,2), 1)/norm(dx0(:,1)) / norm(dx0(:,2)));

T = [
    1 cos(gamma)*norm(dx0(:,2))/norm(dx0(:,1))
    0 sin(gamma)*norm(dx0(:,2))/norm(dx0(:,1))
    ];

Tinv = inv(T);

eta0 = T * xi0;

L_eta = subs(L, xi, T \ eta);

dL_eta = [
    diff(L_eta, eta(1), 1)
    diff(L_eta, eta(2), 1)
    ];
ddL_eta = [
    diff(dL_eta(1,:), eta(1), 1)
    diff(dL_eta(2,:), eta(1), 1)
    diff(dL_eta(2,:), eta(2), 1)
    ];

dx_eta = X * subs(dL_eta, eta, eta0).';
ddx_eta = X * subs(ddL_eta, eta, eta0).';
J0vec = cross(dx_eta(:,1), dx_eta(:,2));
J0 = norm(J0vec);

N0 = 1;

result = 0;

for n = 1 : 4
    c1 = T * Corners(:,n);
    c2 = T * Corners(:,mod(n+1-1,4)+1);
    
    l = c2-c1;
    l = l / norm(l);
    
    d2 = c2-eta0;
    d1 = c1-eta0;
    
    d0 = d2 - l*dot(d2,l);
    d = norm(d0);
    
    t0 = atan2(d0(2), d0(1));
    t1 = atan2(d1(2), d1(1));
    t2 = atan2(d2(2), d2(1));
    if (t1 > t2)
        t2 = t2 + 2*pi;
    end
    
    [theta, w] = gaussquad(nGauss, t1, t2);
    
    Avec = dx_eta * [
        cos(theta).'
        sin(theta)'
        ];
    Bvec = ddx_eta * [
        cos(theta).^2.' / 2
        sin(theta).' .* cos(theta).'
        sin(theta).^2.' / 2
        ];
    J1vec = (cross(ddx_eta(:,1), dx_eta(:,2)) + cross(dx_eta(:,1), ddx_eta(:,2))) * cos(theta).' +...
        (cross(ddx_eta(:,2), dx_eta(:,2)) + cross(dx_eta(:,1), ddx_eta(:,3))) * sin(theta).';
    A = sqrt(dot(Avec, Avec, 1)).';
    
    F2 = J0 * N0 ./ (4*pi*A.^3);
    
    rho_lim = d ./ cos(theta-t0);
    
    %     disp(w);
    
    hold on;
    plot(theta, (-F2 ./ rho_lim));
    
    result = result + w.' * (-F2 ./ rho_lim);
end

fprintf(1, 'result: %.12g\n', result);