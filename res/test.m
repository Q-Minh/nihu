clear;

[mu_lim_c, Teta_c] = quad_match('corner');
[mu_lim_e, Teta_e] = quad_match('edge');
[mu_lim_f, Teta_f] = quad_match('face');

for n = 1 : 10
    disp(n);
    
    [Xi, W] = duffy_general(mu_lim_c, Teta_c, n);
    I_c(n) = integrate(Xi, W, Teta_c);
    N_c(n) = length(W);
    
    [Xi, W] = duffy_general(mu_lim_e, Teta_e, n);
    I_e(n) = integrate(Xi, W, Teta_e);
    N_e(n) = length(W);
    
    [Xi, W] = duffy_general(mu_lim_f, Teta_f, n);
    I_f(n) = integrate(Xi, W, Teta_f);
    N_f(n) = length(W);
end

figure;
hold on;
plot(N_c, log10(abs(I_c/I_c(end)-1)), '*r-');
plot(N_e, log10(abs(I_e/I_e(end)-1)), '*b-');
plot(N_f, log10(abs(I_f/I_f(end)-1)), '*k-');
legend('corner', 'edge', 'facial');
set(gca, 'xscale', 'log');
grid;

% figure;
% plot(Xi(:,1), Xi(:,2), 'b.');
% hold on;
% plot(Xi(:,3), -Xi(:,4), 'r.');
% 
