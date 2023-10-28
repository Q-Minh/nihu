clear;

corners = [0 0 0; 1 0 0; 0 1 0]*.1;
%corners = [0 0 0; .5 0 0; 1 0 0; .5 .5 0; 0 1 0; 0 .5 0];

%x0 = mean(corners, 1);
x0 = [1/3 1/3 0]*.1;

% Distance from each node
R = bsxfun(@minus, corners, x0);
r = sqrt(dot(R, R, 2));
R = bsxfun(@times, R, 1./r);

C = corners - circshift(corners, -1, 1);
c = sqrt(dot(C, C, 2));
C = bsxfun(@times, C, 1./c);

theta = acos(dot(R, circshift(R, -1, 1), 2));
alpha = acos(dot(R, C, 2));

% int(1/r)
I_m1 = r .* sin(alpha) .* log(tan((alpha + theta)/2) ./ tan(alpha/2));
fprintf('int(1/r) = %f\n', sum(I_m1));
%I_stat += r[i] * std::sin(alpha[i]) *
%			std::log(std::tan((alpha[i] + theta[i]) / 2.0) / std::tan(alpha[i] / 2.0)

nC = size(corners, 1);
I_m1_stat_quadr = zeros(nC, 1);
I_m1_dyn_quadr = zeros(nC, 1);
I_m1_full_quadr = zeros(nC, 1);

f = 38;
c = 340;
k = 2*pi*f / c;

K_m1_stat = @(r)(1./r);
K_m1_dyn = @(r)(exp(-1i*k*r/2) * (-1i*k) .* sinc(k*r/2/pi));
K_m1_full = @(r)(exp(-1i*k*r)./r);

for iC = 1 : nC
    mesh.Nodes = [[1; 2; 3],[x0; corners([iC, mod(iC, nC)+1], :)]];
    mesh.Elements = [1 ShapeSet.LinearQuad.Id 1 1 1 1 2 3];
    [mesh.Materials, mesh.Properties] = default_mat_prop();
    [xg, ~, wg] = geo2gauss(mesh, 11);
    dx = bsxfun(@minus, xg, x0);
    dr = sqrt(dot(dx, dx, 2));
    I_m1_stat_quadr(iC) = wg.' * K_m1_stat(dr);
    I_m1_dyn_quadr(iC) = wg.' * K_m1_dyn(dr);
    I_m1_full_quadr(iC) = wg.' * K_m1_full(dr);
end

% int(1/r^3)
I_m3 = (cos(alpha + theta) - cos(alpha)) ./ (r .* sin(alpha));
fprintf('int(1/r^3) = %f\n', sum(I_m3));
%(std::cos(alpha[i] + theta[i]) - std::cos(alpha[i])) / (r[i] * std::sin(alpha[i]))