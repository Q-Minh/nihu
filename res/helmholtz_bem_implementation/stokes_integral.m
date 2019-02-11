function I = stokes_integral(coords, lset, x0, nx0, grad_kernel)

domain = lset.Domain;
xi_nodes = domain.CornerNodes;
nCorners = size(xi_nodes,1);

nG = 30;
[eta, w] = gaussquad(nG);

I = 0;

for n = 1 : nCorners
    xi = (1-eta)/2 * xi_nodes(n,:)  +...
        (1+eta)/2 * xi_nodes(mod(n,nCorners)+1,:);
    dxieta = (xi_nodes(mod(n,nCorners)+1,:) - xi_nodes(n,:)) / 2;
    
    [L, dL] = lset.eval(xi);
    y = L * coords;
    dyeta = dL(:,:,1) * coords * dxieta(1) + ...
        dL(:,:,2) * coords * dxieta(2);
    dl = bsxfun(@times, dyeta, w);
    
    gradG = grad_kernel(x0, y);
    
    I = I + sum(cross(gradG, dl, 2), 1) * nx0.';
end

end
