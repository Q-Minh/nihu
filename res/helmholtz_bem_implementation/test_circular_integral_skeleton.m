corners = [
    0 0 0
    1 0 0
    1 1 0
    0 1 0
    ];

nx0 = [1 0 0];
order = 20;

xvec = -.1 : 1e-2 : 1.1;
xvec = sort(unique([xvec, ...
    -1e-1 : 1e-3 : 1e-1, (-1e-1 : 1e-3 : 1e-1) + 1 ...
    -1e-2 : 1e-4 : 1e-2, (-1e-2 : 1e-4 : 1e-2) + 1 ]));
zvec = [2e-2 : -1e-3 : 0];
g = zeros(length(xvec), length(zvec));

for j = 1 : length(zvec)
    progbar(1, length(zvec), j);
    for i = 1 : length(xvec)
        x0 = [xvec(i) .5 zvec(j)];
        g(i,j) = laplace_3d_hsp_nearly_singular(x0, nx0, corners, order);
    end
end

%%
figure;
plot(xvec, g);
legend(num2str(zvec'));
grid;