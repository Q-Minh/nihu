clear;
clc;

x = sym('x', 'real');
D = sym('D', 'positive');
k = sym('k', 'positive');
n = sym('n', 'positive');
evalin(symengine,'assume(n,Type::Integer)');

g = (k*D+x)^(-n);
N = 13;
kvec = 1 : .1 : 25;

Taylor = sym('Taylor', [1, N+1]);
for i = 0 : N
    Taylor(i+1) = subs(diff(g, x, i), x, 0)/factorial(i);
end

I = int(Taylor .* x.^(0:N), x, -D/2, D/2);
Ianal = int(g, x, -D/2, D/2);

Error = zeros(length(kvec),N+1,3);
for nval = 1 : 3
    In = subs(I,n,nval);
    Inanal = subs(Ianal, n, nval);
    for iK = 1 : length(kvec)
        disp(kvec(iK));
        Error(iK,:,nval) = double(subs(In,k,kvec(iK))/subs(Inanal, k, kvec(iK)));
    end
end

o2 = zeros(length(kvec), 3);
o3 = zeros(length(kvec), 3);
for nval = 1 : 3
    o2(:,nval) = squeeze(sum(abs(cumsum(Error(:,:,nval), 2)-1) > 1e-2,2));
    o3(:,nval) = squeeze(sum(abs(cumsum(Error(:,:,nval), 2)-1) > 1e-3,2));
end

figure; plot(kvec, o2);
figure; plot(kvec, o3);
