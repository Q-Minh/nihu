clear;

%%

XX = [0 0];
YY = [2 1];
Y = YY + [-1 1]/2;
y = Y + [1 -1]/4;
X = XX + [1 1]/2;
x = X + [1 1]/4;

% x = [10-0 0];
% y = [10-10 0];
% X = [10-1.25 -1.25];
% Y = [10-8.75 -1.25];

figure;
plot(x(1), x(2), 'r.', ...
    y(1), y(2), 'b.', ...
    X(1), X(2), 'r*', ...
    Y(1), Y(2), 'b*', ...
    XX(1), XX(2), 'ro', ...
    YY(1), YY(2), 'ko'...
    );
text(x(1), x(2), 'x');
text(y(1), y(2), 'y');
text(X(1), X(2), 'X');
text(Y(1), Y(2), 'Y');
text(XX(1), XX(2), 'XX');
text(YY(1), YY(2), 'YY');
axis equal;

N = 20;
mu = (-N : N)';

k = 1/30;
rxy = x - y;
[~, rxy] = cart2pol(rxy(1), rxy(2));
phi0 = besselh(0, 2, k*rxy) * 2 * -1i/4;

%% analytical functions
rfun = @(r)sqrt(sum(r.*r));
thfun = @(r)atan2(r(2), r(1));

P2Mfun = @(mu, Y, y, k)besselj(mu, k*rfun(Y-y)) .* exp(-1i*mu*thfun(Y-y));
M2Lfun = @(mu, X, Y, k)besselh(mu + mu', 2, k*rfun(X-Y)) .* exp(1i*(mu + mu') * (pi + thfun(X-Y)));
L2Pfun = @(mu, x, X, k)besselj(mu', k*rfun(x-X)) .* exp(-1i*mu'*thfun(x-X));

M2Mfun = @(mu, Psi, Y, k)besselj(mu-mu', k*rfun(Psi-Y)) .* exp(-1i*(mu-mu')*(thfun(Psi-Y)));
L2Lfun = @(mu, X, Xi, k)besselj(mu - mu', k*rfun(X-Xi)) .* exp(1i*(mu - mu')*(pi+thfun(X-Xi)));

%%  matrices
p2m = P2Mfun(mu, Y, y, k);
m2m = M2Mfun(mu, YY, Y, k); 
m2l = M2Lfun(mu, XX, YY, k);
l2l = L2Lfun(mu, X, XX, k);
l2p = L2Pfun(mu, x, X, k);

%%

N = p2m * 2 * -1i/4;
MM = m2m * N;
LL = m2l * MM;
L = l2l * LL;
phi = l2p * L;

err = log10(abs(phi./phi0-1));
if (err > -8)
    error('Baj van');
else
    fprintf(1, 'P2M -> M2L -> L2P works with error %g\n', err);
end

% %% testing L2L and M2M
% phi = L2Pfun(mu, x, X, k) *...
%     L2Lfun(mu, X, Xi, k) *...
%     M2Lfun(mu, Xi, Psi, k) *...
%     M2Mfun(mu, Psi, Y, k) *...
%     P2Mfun(mu, Y, y, k);
% err = log10(abs(phi./phi0-1));
% if (err > -10)
%     error('Baj van');
% else
%     fprintf(1, 'P2M -> M2M -> M2L -> L2L -> L2P works with error %g\n', err);
% end
