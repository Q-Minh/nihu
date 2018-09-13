function [theta, alpha, r] = plane_elem_helper(corners, x0)

R = bsxfun(@minus, corners, x0);
r = sqrt(dot(R, R, 2));
R = bsxfun(@times, R, 1./r);

C = corners - circshift(corners, -1, 1);
c = sqrt(dot(C, C, 2));
C = bsxfun(@times, C, 1./c);

theta = acos(dot(R, circshift(R, -1, 1), 2));
alpha = acos(dot(R, C, 2));

end
