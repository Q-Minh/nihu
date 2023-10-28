function [theta, alpha, r] = plane_elem_helper(corners, x0)
% PLANE_ELEM_HELPER Compute helper quantities on a plane element
%   [THETA, ALPHA, R] = PLANE_ELEM_HELPER(CORNERS, X0) Computes the helpers
%       THETA, ALPHA, R for the plane element defined by CORNERS and the
%       singular point X0. 

R = bsxfun(@minus, corners, x0);
r = sqrt(dot(R, R, 2));
R = bsxfun(@times, R, 1./r);

C = corners - circshift(corners, -1, 1);
c = sqrt(dot(C, C, 2));
C = bsxfun(@times, C, 1./c);

theta = acos(dot(R, circshift(R, -1, 1), 2));
alpha = acos(dot(R, C, 2));

end % of plane_elem_helper
