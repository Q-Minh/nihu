function Z = quad_intersect(X, Y)
%QUAD_INTERSECT compute intersection of two rectangles
% Z = QUAD_INTERSECT computes the intersection domain ot two rectangles.
% Both input rectangles are defined in the 4x2 matrices X and Y. The
% resulting rectangle Z has two columns.

mx1 = min(max(X(:,1)), max(Y(:,1)));
mn1 = max(min(X(:,1)), min(Y(:,1)));
mx2 = min(max(X(:,2)), max(Y(:,2)));
mn2 = max(min(X(:,2)), min(Y(:,2)));

Z = [
    mn1 mn2
    mx1 mn2
    mx1 mx2
    mn1 mx2
    ];

end
