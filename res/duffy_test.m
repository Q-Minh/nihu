% figure;
% hold on;
% plot3(Xi(:,1), Xi(:,2), W, '.b');
% plot3(Xi(:,3), -Xi(:,4), W, '.r');
% % axis equal;
% 

for n = 1 : 15
[Xi, W] = duffy_tria(n, 'face');


r = sqrt((Xi(:,3)-Xi(:,1)).^2 + (Xi(:,4)-Xi(:,2)).^2);
f = 1./ r;
I(n) = W' * f;
end

plot(log10(abs(I / I(15)-1)), 'r*-');

