corners = [
    0 0 0
    1 0 0
    0 .7 0
    ];
xvec = -1 : 1e-2 : 2;
n = [1 1 1];
n = n / norm(n);

order = 40;
for i = 1 : length(xvec)
    x0 = [xvec(i), .3, .1];
    
    g(i) = laplace_3d_slp_nearly_singular(x0, [], corners, order);
    h(i) = laplace_3d_dlp_nearly_singular(x0, [], corners, order);
    ht(i) = laplace_3d_dlpt_nearly_singular(x0, n, corners, order);
    d(i) = laplace_3d_hsp_nearly_singular(x0, n, corners, order);
end

figure;
plot(xvec, g, ...
    xvec, h, ...
    xvec, ht, ...
    xvec, d);