function I = integrate(XI, W, Teta)

mu = XI(:,1:2) + XI(:,3:4)*Teta;
r = sqrt(dot(mu,mu,2));
I = W' * (1./r);

end
