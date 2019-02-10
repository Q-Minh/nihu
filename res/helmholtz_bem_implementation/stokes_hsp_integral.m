function  I = stokes_hsp_integral(corners, lset, x, nx,...
    y0, N0, gradN0, ...
    slp_kernel, grad_kernel)

domain_corners = lset.Domain.CornerNodes;
nVert = size(domain_corners, 1);

I = zeros(size(gradN0,1),1);

% perform contour integrals
order = 40;
nG = ceil((order+1) / 2);  % num of Gaussian points
[eta, w] = gaussquad(nG);       % line quadrature

for i = 1 : nVert
    % line limits in intrinsic coordinates
    xi1 = domain_corners(i,:);
    xi2 = domain_corners(mod(i,nVert)+1,:);
    dxi_eta = (xi2 - xi1) / 2;
    
    xi = (1-eta)/2 * xi1 + (1+eta)/2 * xi2;
    [L, dL] = lset.eval(xi);
    y = L * corners;
    dy = dL(:,:,1) * corners * dxi_eta(1) + ...
        dL(:,:,2) * corners * dxi_eta(2);
    
    G = slp_kernel(x, y);
    gradG = grad_kernel(x, y);
    
    I = I + N0 * dot(nx, w.' * cross(gradG, dy, 2)) ...
        + gradN0 * (...
        bsxfun(@minus, y, y0).' * diag(w) * cross(gradG, dy, 2) * nx.' ...
        - cross(dy, repmat(nx, nG, 1), 2).' * diag(w) * G);
end

% compute surface integral term (assuming full regularity)
[xi, w] = gaussquad2(order, nVert);
[L, dL] = lset.eval(xi);
y = L* corners;
jvec = cross(dL(:,:,1) * corners, dL(:,:,2) * corners, 2);
jac = sqrt(dot(jvec, jvec, 2));
ny = bsxfun(@times, jvec, 1./jac);
gradG = grad_kernel(x, y);
dGnx = -gradG * nx.';
I = I + gradN0 * ny.' * diag(w .* jac) * dGnx;

figure;
plot3(xi(:,1), xi(:,2), bsxfun(@times, (gradN0 * ny.').', diag(jac) * dGnx), '.');
title('Surface integral in Stokes mehod');

% compute spherical angle numerically
dGny = dot(gradG, ny, 2);
Om = -4*pi* (w .* jac).' * dGny;

I = I - (gradN0 * nx.') * Om/(4*pi);

end % of function
