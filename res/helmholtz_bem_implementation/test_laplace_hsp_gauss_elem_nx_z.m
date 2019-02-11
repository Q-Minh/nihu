clear;

coords = [
    -1 -1
    +1 -1
    +1 +1
    -1 +1
    ];

z = 1+sqrt(3)/3;

xi0 = [-1, -1/sqrt(3)];

[theta_lim, theta_0, ref_distance] = plane_elem_helper_mid(coords, xi0);

eta = xi0(2);
xi = xi0(1);

N0 = [
    (1.0 - sqrt(3.0)*xi) * (1.0 - sqrt(3.0)*eta)
    (1.0 + sqrt(3.0)*xi) * (1.0 - sqrt(3.0)*eta)
    (1.0 - sqrt(3.0)*xi) * (1.0 + sqrt(3.0)*eta)
    (1.0 + sqrt(3.0)*xi) * (1.0 + sqrt(3.0)*eta)
    ] / 4.0;

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
    
    Rlim = ref_distance(i) ./ (cos(th - theta_0(i)));
    rlim = sqrt(z^2 + Rlim.^2);
    
    c = cos(th);
    s = sin(th);
    
    int1 = c .* Rlim.^3 ./ (z* rlim.^3);
    
    int2 = -c * z .* (3*Rlim.^2 + 2*z^2) ./ rlim.^3;
    
    int3 = 3*c*z .* (log(rlim + Rlim) - Rlim .* (3*rlim.^2 + Rlim.^2) ./ 3 ./ rlim.^3);
    
    I = I...
        + N0 * (wth.' * int1) ...
        + dN0(:,1) * (wth.' * (c .* int2)) ...
        + dN0(:,2) * (wth.' * (s .* int2)) ...
        + dNxy * (wth.' * (c.*s .* int3));
    
    I0 = I0 + N0 * (wth.' * int1);
end

z

I/4/pi