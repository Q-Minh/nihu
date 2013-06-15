function X = quad_intersect(mu, eta_lim)

X = bsxfun(@plus, -mu, eta_lim);
X = min(max(X,0),1);

end
