clear;
clc;

K = @(x,y)1./abs(bsxfun(@minus, x(:), y(:).'));

%%
Nx = 30;        % number of points in X interval
X = [1 2]';      % X interval
Ny = 20;        % number of points in Y interval
Y = [3 4]';      % Y interval
n = 4;          % expansion length

% generate points
x = linspace(X(1), X(2), Nx)';
y = linspace(Y(1), Y(2), Ny)';

%%
P2M = chebinterp(n, y, Y);
M2L = K(lintrans(chebroots(n), [-1 1]', X),...
    lintrans(chebroots(n), [-1 1]', Y));
L2P = chebinterp(n, x, X)';

%%
sigma = rand(Ny,1); % excitation points

Far = P2M*sigma;
Near = M2L * Far;
f = L2P * Near;

f0 = K(x, y) * sigma;

figure;
plot(x, f, x, f0);
figure;
plot(x, abs(f-f0)./abs(f0))

