function Z = tri_intersect(mu, eta_lim)

s = -mu + min(eta_lim,[],1);

s1 = s(1);
s2 = s(2);

c2 = max(s2,0);
c3 = min(1, 1+s1);
c1 = c2 + max(s1-s2,0);
c4 = c2 + c3 - c1;

Z = [
    c1 c2
    c3 c2
    c3 c4
    ];

end
