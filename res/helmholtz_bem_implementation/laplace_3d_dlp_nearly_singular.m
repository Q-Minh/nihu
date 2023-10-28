function g = laplace_3d_dlp_nearly_singular(x, ~, corners, order)

if nargin < 4
    order = 20;
end

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

z = x(3) - corners(1,3);

[theta_lim, theta_0, ref_distance] = plane_elem_helper_mid(...
    bsxfun(@times, corners, [1 1 0]), ...
    x .* [1 1 0]);

[xi, w] = gaussquad1(order);
N = ShapeSet.LinearLine.eval(xi);

nC = size(corners, 1);
for i = 1 : nC
    idx = [i, mod(i, nC)+1 ];
    
    th = theta_lim(idx);
    if (abs(th(2) - th(1)) > pi)
        [~, j] = min(th);
        th(j) = th(j) + 2*pi;
    end
    
    theta = N * th;
    w_theta = w * diff(th)/2;
    
    R = ref_distance(i) ./ cos(theta - theta_0(i));
    r = sqrt(R.^2 + z^2);
    
    integrand = 1 - abs(z) ./ r;
    
    g = g + sum(integrand .* w_theta);
end

g = g * sign(z) / (4*pi);

end % of function
