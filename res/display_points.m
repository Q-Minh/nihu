function big_test
N = 4;

types = {'face', 'edge', 'corner'};

for i = 1 : length(types)
    [Xi, W] = duffy_tria_tria(N, types{i});
    f(i) = plot_quad(Xi, W);
end

for i = 2 : length(types)
    [Xi, W] = duffy_tria_quad(N, types{i});
    f(i) = plot_quad(Xi, W);
end

end

function f = plot_quad(Xi, W)
f = figure;
subplot(1,2,1);
plot(Xi(:,1), Xi(:,2), 'b.');
xlabel('x_1');
ylabel('x_2');
axis equal;
title(sprintf('N = %d', length(W)));
subplot(1,2,2);
plot(Xi(:,3), Xi(:,4), 'r.');
xlabel('y_1');
ylabel('y_2');
axis equal;
title(sprintf('Sum(w) = %g', sum(W)));
end
