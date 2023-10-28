function convergence_test

figure;
hold on;
p = zeros(0,1);
N = 15;

data = {
    'face',   @duffy_tria_tria, 23, [0 0;  1 0; 0  1], 'b*-'
    'edge',   @duffy_tria_tria, 23, [0 0;  1 0; 0 -1], 'bo-'
    'corner', @duffy_tria_tria, 23, [0 0; -1 0; 0 -1], 'b.-'
    ...
    'edge',   @duffy_tria_quad, 24, [0 0; 1 0; 1 -1; 0 -1], 'ro-'
    'corner', @duffy_tria_quad, 24, [0 0; -1 0; -1 -1; 0 -1], 'r.-'
    };


for d = 1 : size(data,1)
    I = zeros(N,1);
    num = zeros(N,1);
    for n = 1 : N
        fun = data{d,2};
        [Xi, W] = fun(n, data{d,1});
        C1 = [0 0; 1 0; 0 1];
        C2 = data{d,4};
        [x, wx] = quadrature_map(Xi(:,1:2), sqrt(W), C1, 23);
        [y, wy] = quadrature_map(Xi(:,3:4), sqrt(W), C2, data{d,3});
        num(n) = length(W);
        I(n) = integrate(x, y, abs(wx.*wy));
        p(end+1,1) = plot(num, log10(abs(I/I(N)-1)), data{d,5});
    end
end
set(gca, 'xscale', 'log');
grid;

end

function I = integrate(x, y, w)
rvec = y - x;
r = sqrt(dot(rvec, rvec, 2));
I = w' * (exp(-1i*r)./r.^2);
end
