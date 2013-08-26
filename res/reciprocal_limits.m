clear;
clc;

%% computation parameters
N = 13;                 % maximal polynomial order
kvec = (1 : .1 : 50)';     % relative distance vector
n_max = 3;              % maximal inverse order

%% declare symbolic Green's function
x = sym('x', 'real');
D = sym('D', 'positive');
k = sym('k', 'positive');
n = sym('n', 'positive');
evalin(symengine,'assume(n,Type::Integer)');
g = (k*D+x)^(-n);

%% Taylor expansion
disp('Computing Taylor expansion');
Taylor = sym('Taylor', [1, N+1]);
for i = 0 : N
    Taylor(i+1) = subs(diff(g, x, i), x, 0)/factorial(i);
end

%% Integrate the expansion and the analytical function
disp('Integrating');
I = int(Taylor .* x.^(0:N), x, -D/2, D/2);
Ianal = int(g, x, -D/2, D/2);

%% Compute relative error
fprintf(1, 'Computing relative errors\n');
Error = zeros(length(kvec),N+1,n_max);
nb = 0;
for nval = 1 : n_max
    In = subs(I,n,nval);
    Inanal = subs(Ianal, n, nval);
    for iK = 1 : length(kvec)
        fprintf(1, repmat('\b', 1, nb));
        nb = fprintf(1, 'n = %d - k = %g', nval, kvec(iK));
        Error(iK,:,nval) = double(...
            subs(In,k,kvec(iK))/subs(Inanal, k, kvec(iK)));
    end
end

%% Compute polynomial orders for given accuracy
disp('Computing polynomial orders');
o2 = zeros(length(kvec), n_max);
o3 = zeros(length(kvec), n_max);
o4 = zeros(length(kvec), n_max);
for nval = 1 : n_max
    o2(:,nval) = squeeze(sum(abs(cumsum(Error(:,:,nval), 2)-1) > 1e-2,2));
    o3(:,nval) = squeeze(sum(abs(cumsum(Error(:,:,nval), 2)-1) > 1e-3,2));
    o4(:,nval) = squeeze(sum(abs(cumsum(Error(:,:,nval), 2)-1) > 1e-4,2));
end

%% display results
figure; plot(kvec, o2);
title('accuracy: 1e-2');
xlabel('relative distance [-]');
ylabel('required polynomial order');
legend('1/r', '1/r^2', '1/r^3');

figure; plot(kvec, o3);
title('accuracy: 1e-3');
legend('1/r', '1/r^2', '1/r^3');
xlabel('relative distance [-]');
ylabel('required polynomial order');

figure; plot(kvec, o4);
title('accuracy: 1e-4');
legend('1/r', '1/r^2', '1/r^3');
xlabel('relative distance [-]');
ylabel('required polynomial order');


%%
clc;
res = {o2 o3 o4};
for r = 1 : length(res)
    fprintf(1, 'eps = %g\n', 10^(-r-1));
    for n = 1 : 3
        fprintf(1, '\tn = %g\n', n);
        f = find(diff(res{r}(:,n)));
        disp([kvec(f) res{r}(f,n)]);
    end
end