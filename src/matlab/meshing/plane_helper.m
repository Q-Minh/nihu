function [R, alpha, theta, theta0] = plane_helper(coords, x0)

nc = size(coords,1);

theta = nan(nc,1);
alpha = nan(nc,1);
R = nan(nc,1);
theta0 = nan(nc,1);

% loop through triangles
for ic = 1 : nc
    % vector from center to corner
    V = coords(ic,:) - x0;
    next = ic + 1;
    if next > nc
        next = 1;
    end
    Vnext = coords(next,:) - x0;
    t = coords(next,:) - coords(ic,:);
    
    R(ic) = norm(V);
    alpha(ic) = acos(dot(V, -t, 2) / norm(V) / norm(t));
    theta(ic) = acos(dot(V, Vnext, 2) / norm(V) / norm(Vnext));
    
    theta0(ic) = atan2(V(2), V(1));
end

end % of function
