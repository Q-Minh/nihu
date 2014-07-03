N = 4;
[x,w] = gaussquad(N,0,1);
w = w*w';
w = w(:);
w = w*w';
w = w(:);


[x1,x2,x3,x4] = ndgrid(x,x,x,x);
x1 = x1(:);
x2 = x2(:);
x3 = x3(:);
x4 = x4(:);

J = x3 .* x1.^3;
xi1 = x1;
xi2 = x1.*x2;
eta1 = x1.*x3;
eta2 = x1.*x3.*x4;

Xi1 = [xi1; eta1];
Xi2 = [xi2; eta2];
Eta1 = [eta1; xi1];
Eta2 = [eta2; xi2];
W = [w.*J; w.*J];

figure;
hold on;
plot3(Xi1, Xi2, W, 'b.');
plot3(-Eta1, Eta2, W,'r.');