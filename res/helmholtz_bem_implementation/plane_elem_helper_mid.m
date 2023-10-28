function [theta_lim, theta_0, ref_distance] = plane_elem_helper_mid(corners, x0)
% PLANE_ELEM_HELPER Compute helper quantities on a plane element
%   [THETA, ALPHA, R] = PLANE_ELEM_HELPER(CORNERS, X0) Computes the helpers
%       THETA, ALPHA, R for the plane element defined by CORNERS and the
%       singular point X0.

nCorners = size(corners,1);
theta_lim = nan(nCorners,1);
theta_0 = nan(nCorners,1);
ref_distance = nan(nCorners,1);

for i = 1 : nCorners
    c1 = corners(i,:);
    c2 = corners(mod(i, nCorners)+1,:);
    
    l = c2-c1;
    l = l / norm(l);
    
    d1 = c1 - x0;
    d0 = d1 - l * dot(d1, l, 2);
    
    theta_lim(i) = atan2(d1(2), d1(1));
    theta_0(i) = atan2(d0(2), d0(1));
    ref_distance(i) = norm(d0);
end

end % of function
