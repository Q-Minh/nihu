% clear;

corners = [
    0 0 0
    1 0 0
    0 .7 0
    ];
xvec = -1 : 1e-2 : 2;
n = [1 1 1];
n = n / norm(n);

order = 10;
% for i = 1 : length(xvec)
    x0 = [.3, .3, .3];
    
    k = 1;
    
    g = helmholtz_3d_slp_nearly_singular(x0, [], corners, k, order);
%     h = laplace_3d_dlp_nearly_singular(x0, [], corners, order);
%     ht = laplace_3d_dlpt_nearly_singular(x0, n, corners, order);
%     d = laplace_3d_hsp_nearly_singular(x0, n, corners, order);
% end

% % figure;
% hold on;
% plot(xvec, g, '-', ...
%     xvec, h, '-', ...
%     xvec, ht, '-', ...
%     xvec, d, '-');