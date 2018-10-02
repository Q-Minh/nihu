function [xi, niter] = be_invmap(coords, x0, shapeset)

coords = coords.';

dim = shapeset.Domain.Space.Dimension;

xi0 = zeros(1, dim);
zeta = 0;


% rotates pi/2 in positive direction if applied from the left

% Newton Raphson cycle
eps = 1e-8;
niter = 1;
while true
    [N, dN, ddN] = shapeset.eval(xi0);
    x = N * coords;
    
    if (dim == 1)
        dx = dN * coords;
        ddx = ddN * coords;
        
        T = [0 -1; 1 0];
        n = (T * dx.').';
        
        dn = (T * ddx.').';
        
    elseif (dim == 2)
        dxxi = dN(:,:,1) * coords;
        dxeta = dN(:,:,2) * coords;
        n = cross(dxxi, dxeta, 2);
        
        dxxixi = ddN(:,:,1,1) * coords;
        dxxieta = ddN(:,:,1,2) * coords;
        dxetaeta = ddN(:,:,2,2) * coords;
        
        dnxi = cross(dxxixi, dxeta, 2) + cross(dxxi, dxxieta, 2);
        dneta = cross(dxxieta, dxeta, 2) + cross(dxxi, dxetaeta, 2);
    end
    
    f = (x + n*zeta) - x0;
    if norm(f) < eps
        break;
    end
    
    if (dim == 1)
    df = [
        dx + dn*zeta
        n
        ];
    elseif (dim == 2)
    df = [
        dxxi + dnxi*zeta
        dxeta + dneta*zeta
        n
        ];
    end
    xi = [xi0 zeta] - (df.' \ f.').';
    
    xi0 = xi(1:end-1);
    zeta = xi(end);

    niter = niter + 1;
    if (niter > 100)
        break;
    end

end

xi = [xi0 zeta];

end
