function Phi = iePhij3(P,alpha,beta)

z = linspace(-1,1-2/P,P);
Phi = zeros(P);
for c1 = 1:P
    for c2 = 1:P
       s = z(c1);
       Phi(c2,c1) = 1/2*(1-s)*jacpoly(c2-1,alpha,beta,s);
    end
end
Phi = Phi.';