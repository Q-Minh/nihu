syms z
N = 14;
nu = 0;
Pnu = P(nu,z,N);
Qnu = Q(nu,z,N);

U = sqrt(Pnu^2 + Qnu^2);
Phi = atan(Qnu/Pnu);

syms x;
U = subs(U, z, 1/x);
Phi = subs(Phi, z, 1/x);

U = taylor(U, x, 'Order', N);
Phi = taylor(Phi, x, 'Order', N);

for n = 0 : N
    Unum(n+1) = limit(diff(U,x,n)/factorial(n), x, 0);
    Phinum(n+1) = limit(diff(Phi,x,n)/factorial(n), x, 0);
end
