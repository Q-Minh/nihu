function I = integrate(XI, W)

mu = XI(:,3:4) - XI(:,1:2);
r = sqrt(dot(mu,mu,2));
I = W' * (1./r);

end
