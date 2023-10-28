corners = [
    0 0 0
    1 0 0
    0 1 0
    ];

order = 10;

x = bsxfun(@plus, [.3 .3 0], (-1e-3 : 1e-5 : 1e-3)'/1e3 * [0 0 1]);
% x = x + 1e-8 * randn(size(x));
N = size(x, 1);
nx = [1 1 1];
nx = nx / norm(nx);

%%
for i = 1 : N
    slp(i) = laplace_3d_slp_nearly_singular(x(i,:), [], corners, order);
    dlp(i) = laplace_3d_dlp_nearly_singular(x(i,:), [], corners, order);
    dlpt(i) = laplace_3d_dlpt_nearly_singular(x(i,:), nx, corners, order);
    hsp(i) = laplace_3d_hsp_nearly_singular(x(i,:), nx, corners, order);
end

%%
figure;
subplot(2,2,1);
plot(x(:,3), slp);
subplot(2,2,2);
plot(x(:,3), dlp);
subplot(2,2,3);
plot(x(:,3), dlpt);
subplot(2,2,4);
plot(x(:,3), hsp);


