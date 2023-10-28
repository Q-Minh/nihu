function N = quad1_gauss_shapeset(xivec)

xi = xivec(:,1);
eta = xivec(:,2);

N = [
    (1.0 - sqrt(3.0)*xi) .* (1.0 - sqrt(3.0)*eta), ...
    (1.0 + sqrt(3.0)*xi) .* (1.0 - sqrt(3.0)*eta), ...
    (1.0 - sqrt(3.0)*xi) .* (1.0 + sqrt(3.0)*eta), ...
    (1.0 + sqrt(3.0)*xi) .* (1.0 + sqrt(3.0)*eta)
    ] / 4;

end
