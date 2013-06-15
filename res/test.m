clear;

[mu_lim_c, eta_lim_c] = quad_match('corner');
[mu_lim_e, eta_lim_e] = quad_match('edge');
[mu_lim_f, eta_lim_f] = quad_match('face');

for n = 1 : 6
    disp(n);
    
    [Xi, W] = duffy_general(mu_lim_f, eta_lim_f, n, @quad_intersect, @quad_map);
    I_f(n) = integrate(Xi, W);
    N_f(n) = length(W);

    [Xi, W] = duffy_general(mu_lim_e, eta_lim_e, n, @quad_intersect, @quad_map);
    I_e(n) = integrate(Xi, W);
    N_e(n) = length(W);
    
    [Xi, W] = duffy_general(mu_lim_c, eta_lim_c, n, @quad_intersect, @quad_map);
    I_c(n) = integrate(Xi, W);
    N_c(n) = length(W);
end

figure;
hold on;
plot(N_c, log10(abs(I_c/I_c(end)-1)), '*r-');
plot(N_e, log10(abs(I_e/I_e(end)-1)), '*b-');
plot(N_f, log10(abs(I_f/I_f(end)-1)), '*k-');
legend('corner', 'edge', 'facial');
set(gca, 'xscale', 'log');
grid;

figure;
plot(Xi(:,1), Xi(:,2), 'b.');
hold on;
plot(Xi(:,3), Xi(:,4), 'r.');

