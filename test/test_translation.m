clear;
%%
p = 20;
[S, W, perm, model] = spherequad(p);

%%
x = -100:.1:100;
k = 100;
D = [0 5 0];
L = round(3*log10(k*norm(D + [max(x) 0 0])));


%%
M = zeros(size(S,2), length(x));
for i = 1 : length(x)
    disp(i);
M(:,i) = translation(k, D + [x(i) 0 0], L, S);
end

%%
figure;
surf(20*log10(abs(M))); shading interp;


%%

plot_mesh(model, M(:,end));
shading interp;