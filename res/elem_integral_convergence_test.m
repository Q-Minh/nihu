syms r k
Gs = exp(-1i*k*r)/r;
G = matlabFunction(Gs);
dG = matlabFunction( diff(Gs, r));
ddG = matlabFunction(diff(diff(Gs, r), r));

clear r k;
clc;

nodes = [
    0 0 0
    .1 0 0
    0 .1 0
    0 -.1 0
    ];
elems = [
    1 2 3
    1 2 4
    ];

% collocation point
x0 = mean(nodes(elems(1,:),:));
nx = [0 0 -1];

f = 10;
c = 340;
k = 2*pi*f/c;

ordervec4 = 1 : 40;
for nO = 1 : length(ordervec4)
    order = ordervec4(nO);
    
    [xi, w] = gaussquad2(order, 4);
    
    [N, dN] = ShapeSet.LinearQuad.eval(xi);
    elem = elems(2,[1 2 3 3]);
    xg = N * nodes(elem,:);
    xgxi = dN(:,:,1) * nodes(elem,:);
    xgeta = dN(:,:,2) * nodes(elem,:);
    jvec = cross(xgxi, xgeta, 2);
    jac = sqrt(dot(jvec, jvec, 2));
    ny = bsxfun(@times, jvec, 1./jac);
    
    rvec = bsxfun(@minus, xg, x0);
    r = sqrt(dot(rvec, rvec, 2));
    
    rdnx = bsxfun(@times, -rvec * nx.', 1./r);
    rdny = bsxfun(@times, dot(rvec, ny, 2), 1./r);
    
    kernel = (ddG(k,r) - dG(k, r)./r) .* rdnx .* rdny - dG(k,r)./r .* (ny * nx.');
    I4(nO) = w' * diag(jac) *  kernel;
end


ordervec3 = 1 : 9;
for nO = 1 : length(ordervec3)
    order = ordervec3(nO);
    
    [xi, w] = gaussquad2(order, 3);
    
    [N, dN] = ShapeSet.LinearTria.eval(xi);
    elem = elems(2,[1 2 3]);
    xg = N * nodes(elem,:);
    xgxi = dN(:,:,1) * nodes(elem,:);
    xgeta = dN(:,:,2) * nodes(elem,:);
    jvec = cross(xgxi, xgeta, 2);
    jac = sqrt(dot(jvec, jvec, 2));
    ny = bsxfun(@times, jvec, 1./jac);
    
    rvec = bsxfun(@minus, xg, x0);
    r = sqrt(dot(rvec, rvec, 2));
    
    rdnx = bsxfun(@times, -rvec * nx.', 1./r);
    rdny = bsxfun(@times, dot(rvec, ny, 2), 1./r);
    
    kernel = (ddG(k,r) - dG(k, r)./r) .* rdnx .* rdny - dG(k,r)./r .* (ny * nx.');
    I3(nO) = w' * diag(jac) *  kernel;
end

figure;
plot(ordervec3, log10(abs((I3-I4(end))/I4(end))), ...
    ordervec4, log10(abs((I4-I4(end))/I4(end))));
