nExp = 5;

%% 1D
I = [-1; 1];
x0 = chebroots(nExp, 1);
x = (rand(100,1)-.5)*2;
S = chebinterp(nExp, x, I);
z0 = dot(x0,x0,2);
z = S' * z0;

figure;
plot(x0(:,1), z0, 'b.', x(:,1), z, 'r.');

%% 2D
I = [-1 -1; 1 1];
x0 = chebroots(nExp, 2);
x = (rand(100,2)-.5)*2;
S = chebinterp(nExp, x, I);
z0 = dot(x0,x0,2);
z = S' * z0;

figure;
plot3(x0(:,1), x0(:,2), z0, 'b.', x(:,1), x(:,2), z, 'r.');
