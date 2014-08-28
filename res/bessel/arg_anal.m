syms z
N = 14;
nu = 0;
Pnu = P(nu,z,N);
Qun = Q(nu,z,N);

U = sqrt(Pnu^2 + Qun^2);
Phi = atan(Qun/Pnu);

syms x;
U = subs(U, z, 1/x);
Phi = subs(Phi, z, 1/x);

U = taylor(U, x, 'Order', N);
Phi = taylor(Phi, x, 'Order', N);

for n = 1 : N+1
    Unum(n+1) = limit(diff(U,x,n)/factorial(n), x, 0);
    Phinum(n+1) = limit(diff(Phi,x,n)/factorial(n), x, 0);
end
