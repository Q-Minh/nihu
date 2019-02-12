clear;

coords = [
    -1 -1
    +1 -1
    +1 +1
    -1 +1
    ];

xi0 = [-1/sqrt(3), -1/sqrt(3)];

[theta_lim, theta_0, ref_distance] = plane_elem_helper_mid(coords, xi0);

N0 = [1 0 0 0].';

eta = xi0(2);
xi = xi0(1);
dN0 = 	[
    -sqrt(3.0) * (1.0 - sqrt(3.0)*eta) , (1.0 - sqrt(3.0)*xi) * -sqrt(3.0)
    +sqrt(3.0) * (1.0 - sqrt(3.0)*eta) , (1.0 + sqrt(3.0)*xi) * -sqrt(3.0)
    -sqrt(3.0) * (1.0 + sqrt(3.0)*eta) , (1.0 - sqrt(3.0)*xi) * +sqrt(3.0)
    +sqrt(3.0) * (1.0 + sqrt(3.0)*eta) , (1.0 + sqrt(3.0)*xi) * +sqrt(3.0)
    ] / 4.0;

dNxy = 	[
    -sqrt(3.0) * (- sqrt(3.0))
    +sqrt(3.0) * (- sqrt(3.0))
    -sqrt(3.0) * (+ sqrt(3.0))
    +sqrt(3.0) * (+ sqrt(3.0))
    ] / 4.0;

I = zeros(4,1);
I0 = zeros(4,1);
ngauss = 10;
for i = 1 : 4
    t(1) = theta_lim(i);
    t(2) = theta_lim(mod(i,4)+1);
    
    if (abs(t(2)-t(1)) > pi)
        [~, j] = min(t);
        t(j) = t(j) + 2*pi;
    end
    
    [th, wth] = gaussquad(ngauss, t(1), t(2));
    
    rlim = ref_distance(i) ./ (cos(th - theta_0(i)));
    
    I = I -N0 * (wth.' * (1./rlim)) ...
        + dN0(:,1) * (wth.' * (cos(th) .* log(rlim))) ...
        + dN0(:,2) * (wth.' * (sin(th) .* log(rlim))) ...
        + dNxy * (wth.' * (cos(th).*sin(th) .* rlim));
    
    I0 = I0 -N0 * (wth.' * (1./rlim));
end
