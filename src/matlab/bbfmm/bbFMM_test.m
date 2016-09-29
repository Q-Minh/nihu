clear;
clc;

K = @(x,y)laplace_kernel(x,y);

%%
Nx = 30;        % number of points in X interval
X = [1 2]';      % X interval
Ny = 20;        % number of points in Y interval
Y = [3 4]';      % Y interval
nCheb = 8;          % expansion length

% generate points
x = linspace(X(1), X(2), Nx)';  % receivers
y = linspace(Y(1), Y(2), Ny)';  % sources

%%
P2M = chebinterp(nCheb, y, Y);
xNodes = lintrans(chebroots(nCheb), [-1 1]', X);
yNodes = lintrans(chebroots(nCheb), [-1 1]', Y);
M2L = K(xNodes, yNodes);
L2P = chebinterp(nCheb, x, X)';

%%
sigma = rand(Ny,1); % excitation points

%% fmm
Multi = P2M*sigma;  % anterpolation
Local = M2L * Multi; % kernel evaluation
f = L2P * Local; % interpolation

%% conventional direct summation
f0 = K(x, y) * sigma;

%% plot results and check error
figure;
formatfig;
subplot(2,1,1);
plot(x, f, x, f0);
legend('fmm', 'conventional', 'location', 'northwest');
xlabel('x');
subplot(2,1,2);
plot(x, abs(f-f0)./abs(f0))
xlabel('x');
legend('error');

