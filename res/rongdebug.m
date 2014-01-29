clear;

close all;

order = 9;
nGauss = floor(order/2+1);

xi = sym('xi', [2 1]);
eta = sym('eta', [2 1]);

ReferenceCorners = [
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

L_eta = simple(subs(L, xi, T \ eta));

dL_eta = [
    diff(L_eta, eta(1), 1)
    diff(L_eta, eta(2), 1)
    ];
ddL_eta = [
    diff(dL_eta(1,:), eta(1), 1)
    diff(dL_eta(2,:), eta(1), 1)
    diff(dL_eta(2,:), eta(2), 1)
    ];

dx0_eta = X * subs(dL_eta, eta, eta0).';
ddx0_eta = X * subs(ddL_eta, eta, eta0).';
J0vec_eta = cross(dx0_eta(:,1), dx0_eta(:,2));
J0_eta = norm(J0vec_eta);
n0_eta = J0vec_eta / J0_eta;

N0 = 1;

result = 0;

for n = 1 : 4
    c1 = T * ReferenceCorners(:,n);
    c2 = T * ReferenceCorners(:,mod(n+1-1,4)+1);
    
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
    
    [theta, w] = gaussquad(nGauss, double(t1), double(t2));
    
    Avec = double(dx0_eta) * [
        cos(theta).'
        sin(theta)'
        ];
    A = sqrt(dot(Avec, Avec, 1)).';
    Bvec = double(ddx0_eta) * [
        cos(theta).^2.' / 2
        sin(theta).' .* cos(theta).'
        sin(theta).^2.' / 2
        ];
    J1vec_eta = (cross(double(ddx0_eta(:,1)), double(dx0_eta(:,2))) + cross(double(dx0_eta(:,1)), double(ddx0_eta(:,2)))) * cos(theta).' +...
        (cross(double(ddx0_eta(:,2)), double(dx0_eta(:,2))) + cross(double(dx0_eta(:,1)), double(ddx0_eta(:,3)))) * sin(theta).';
    
    F2 = double(J0_eta * N0 ./ (4*pi*A.^3));
    F1 = double(J1vec_eta.' * n0_eta) ./ (4*pi*A.^3) -...
        double(3*(n0_eta.' * J0vec_eta) * dot(Avec, Bvec, 1).') ./ (4*pi*A.^5);
    
    rho_lim = double(d ./ cos(theta-t0));
    
    %     disp(w);
    
    figure(1);
    hold on;
    plot(theta, (-F2 ./ rho_lim));
    
    figure(2);
    hold on;
    plot(theta, (F1 .* log(rho_lim)));
    
    result = result + w.' * (F1 .* log(rho_lim) -F2 ./ rho_lim);
end

fprintf(1, 'result: %.12g\n', result);
