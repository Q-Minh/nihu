function J = myBesselJ(nu,z)
arg = z-(nu/2+1/4)*pi;
N = 5;
J = sqrt(2./(pi*z)) .* (cos(arg) .* P(nu,z,N) - sin(arg) .* Q(nu,z,N));
end
