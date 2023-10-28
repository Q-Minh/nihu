function g = nearly_singular_laplace_3d_planar_constant(kernel_str,...
    x, nx, corners, order)

g = 0;

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
nx = (T \ nx.').';

z = x(3) - corners(1,3);

% compute helper quantities for angular integration
[theta_lim, theta_0, ref_distance] = plane_elem_helper_mid(...
    bsxfun(@times, corners, [1 1 0]), ...
    x .* [1 1 0]);

[xi, w] = gaussquad1(order);
N = ShapeSet.LinearLine.eval(xi);

thlims = zeros(0,3);

nC = size(corners, 1);
for i = 1 : nC
    idx = [i, mod(i, nC)+1 ];
    
    th = [theta_lim(idx); theta_0(i)];
    while (abs(max(th) - min(th)) > pi)
        [~, j] = min(th);
        th(j) = th(j) + 2*pi;
    end
    
    th0 = th(3);
    th = th(1:2);
    
    if th0 > min(th) && th0 < max(th)
        thlims = [
            thlims
            th(1) th0 i
            th0 th(2) i
            ];
    else
        thlims = [
            thlims
            th(1) th(2) i
            ];
    end
end

for j = 1 : size(thlims,1)
    i = thlims(j,3);
    th = thlims(j,1:2).';
    
    theta = N * th;
    w_theta = w * diff(th)/2;
    
    R = ref_distance(i) ./ cos(theta - theta_0(i));
    r = sqrt(R.^2 + z^2);
    
    c = cos(theta);
    s = sin(theta);
    
    rvec = [R.*c R.*s];
    rvec(:,3) = -z;
    
    switch lower(kernel_str)
        case 'slp'
            integrand = r - abs(z);
        case 'dlp'
            integrand = 1 - abs(z) ./ r;
        case 'dlpt'
            integrand = (-rvec * nx.') ./ r ...
                + (nx(1) * c + nx(2) * s) .* log( (R + r) ) ...
                - nx(3) * sign(z);
        case 'hsp'
            rdnx = -(rvec * nx.') ./ r;
            if abs(z) > 1e-12
                integrand = -R.^2 ./ r.^2 ./ z .* rdnx;
            else
                integrand = -nx(3)./R;
            end
        otherwise
            error('invalid kernel string: %s\n', kenrel_str);
    end
    
    g = g + sum(integrand .* w_theta);
end

g = g / (4*pi);

end % of function
