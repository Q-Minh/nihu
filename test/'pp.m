load data/Time

%% plot used time as function of DOF
f(1) = figure;
plot(N, Time(:,:,1));
xlabel('DOF [-]');
ylabel('t_{tree}');
legend({'full', 'symm I', 'symm II'}, 'Location', 'NorthWest');
formatfig;

f(2) = figure;
plot(N, Time(:,:,2));
xlabel('DOF [-]');
ylabel('t_{init}');
legend({'full', 'symm I', 'symm II'}, 'Location', 'NorthWest');
formatfig;

f(3) = figure;
plot(N, Time(:,:,3));
xlabel('DOF [-]');
ylabel('t_{iter}');
legend({'full', 'symm I', 'symm II'}, 'Location', 'NorthWest');
formatfig;

same_size(f);
font_size(f, 12);
linewidth(f, 1);

figure(1);
print -depsc t_tree.eps
figure(2);
print -depsc t_init.eps
figure(3);
print -depsc t_iter.eps