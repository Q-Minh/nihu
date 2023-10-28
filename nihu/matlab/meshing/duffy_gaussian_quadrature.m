function [xi, w] = duffy_gaussian_quadrature(domain, order, idx)
switch domain
    case Domain.Tria
        [xi, w] = gaussian_quadrature(Domain.Quad, order);
        x0 = Domain.Tria.CornerNodes;
        x = x0(mod(idx+[0 0 1 2]-1, 3) + 1, :);
        [L, dL] = ShapeSet.LinearQuad.eval(xi);
        xi = L * x;
        dxxi = dL(:,:,1) * x;
        dxeta = dL(:,:,2) * x;
        j = dxxi(:,1) .* dxeta(:,2) - dxxi(:,2) .* dxeta(:,1);
        w = w .* j;
end
end
