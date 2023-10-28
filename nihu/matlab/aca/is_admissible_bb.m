function b = is_admissible_bb(C1, C2, eta)
bb1 = C1.bb;
bb2 = C2.bb;
d1 = norm(diff(bb1));
d2 = norm(diff(bb2));
dist = sqrt(...
    min(min(abs(bsxfun(@minus, bb1(:,1), bb2(:,1)'))))^2 +...
    min(min(abs(bsxfun(@minus, bb1(:,2), bb2(:,2)'))))^2);
b = min(d1, d2) < eta * dist;
end
