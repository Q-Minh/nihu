clear;
clc;

for n = 1:6
    disp(n);
    
    [Xi, W] = duffy_general(n, 'face', @tri_match, @tri_intersect, @tri_map);
    I_f(n) = integrate(Xi, W);
    N_f(n) = length(W);
    
    [Xi, W] = duffy_general(n, 'corner', @tri_match, @tri_intersect, @tri_map);
    I_c(n) = integrate(Xi, W);
    N_c(n) = length(W);
end

figure;
hold on;
plot(N_c, log10(abs(I_c/I_c(end)-1)), '*r-');
plot(N_f, log10(abs(I_f/I_f(end)-1)), '*k-');
legend('corner', 'facial');
set(gca, 'xscale', 'log');
grid;

% figure;
% plot(Xi(:,1), Xi(:,2), 'b.');
% hold on;
% plot(Xi(:,3), Xi(:,4), 'r.');
