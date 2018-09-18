function g = laplace_3d_dlp_nearly_singular(x, ~, corners, order)

if nargin < 4
    order = 20;
end

% transform into local system
T(:,1) = corners(2,:) - corners(1,:);
T(:,2) = corners(3,:) - corners(1,:);
T(:,3) = cross(T(:,1), T(:,2));
T(:,2) = cross(T(:,3), T(:,1));
for i = 1 : 3
    T(:,i) = T(:,i) / norm(T(:,i));
end

corners = (T \ corners.').';
x = (T \ x.').';

[xi, w] = gaussquad1(order);
[N, dN] = ShapeSet.LinearLine.eval(xi);

g = 0;
for i = 1 : 3
    idx = [i, mod(i, 3)+1 ];
    
    y = N * corners(idx,:);
    dyxi = dN * corners(idx,:);
    
    rvec = bsxfun(@minus, y, x);
    r = sqrt(dot(rvec, rvec, 2));
    
    z = -rvec(:,3);
    
    Rvec = rvec(:,1:2);
    R = sqrt(dot(Rvec, Rvec, 2));
    
    integrand = z./abs(z) -z./r;
    
    dtheta = (Rvec(:,1) .* dyxi(:,2) - Rvec(:,2) .* dyxi(:,1)) ./ R.^2;
    
    g = g + sum(integrand .* dtheta .* w);
end

g = g / (4*pi);

end
