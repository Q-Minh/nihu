%% derive p2m and m2p
clear;
syms z zi
K = log(z-zi);
L = 10;
for k = 0 : L
   M2P(k+1,1) = subs(diff(K, zi, k), zi, 0) / factorial(k) * zi^k;
end
m2p(1,1) = log(z);
for k = 1 : L
    m2p(k+1,1) = 1/z^k;
end
p2m = M2P ./ m2p;

%% derive m2m from p2m
clear;
syms zi z0;
M2P(1,1) = sym(1);
L = 10;
for l = 1 : L
    M2P(l+1,1) = -zi^l/l;
end
for k = 0 : L
    G(:,k+1) = subs(diff(M2P, zi, k), zi, z0) * (zi-z0)^k/factorial(k);
end
m2p(1,1) = sym(1);
for k = 1 : L
    m2p(1,k+1) = -(zi-z0)^k/k;
end
for l = 0 : L
    G(l+1,:) = G(l+1,:) ./ m2p;
end

%%
clear;
syms z z0
L = 10;
for k = 0 : L
    L2P(1,k+1) = z^k;
end
for l = 0 : L
    G(l+1,:) = subs(diff(L2P, z, l), z, z0) * (z-z0)^l / factorial(l);
end
for l = 0 : L
    l2p(l+1,1) = (z-z0)^l;
end
for k = 0 : L
    G(:,k+1) = G(:,k+1) ./ l2p;
end
