function [xi, err, niter] = be_invmap(coords, x0, shapeset, tol, maxiter)
%be_invmap inverse mapping on surface elements
%   [xi, err, niter] = be_invmap(coords, x0, shapeset, tol, maxiter)
%   computes the intrinsic coordinate of point x0 over the element defined
%   by coordinates coord and the geometrical shape set shapeset.
%   xi(1:end-1) contains the intrinsic coordinates, and xi(end) contains
%   the normal distance from the element. The algorithm uses a Newton
%   Raphson iteration. The number of iterations is returned in niter.

% default max number of iterations
if nargin < 5
    maxiter = 100;
end
% default tolerance
if nargin < 4
    tol = 1e-6;
end

coords = coords.';

dim = shapeset.Domain.Space.Dimension;

xi0 = zeros(1, dim);
zeta = 0;

% Newton Raphson cycle
for niter = 1 : maxiter
    [N, dN, ddN] = shapeset.eval(xi0);
    x = N * coords;
    
    if (dim == 1)
        dx = dN * coords;
        ddx = ddN * coords;
        T = [0 -1; 1 0];    % rotates pi/2 in positive direction
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
    err = norm(f);
    if err < tol
        break
    end
    
    if dim == 1
        df = [dx + dn*zeta; n];
    elseif dim == 2
        df = [
            dxxi + dnxi*zeta
            dxeta + dneta*zeta
            n
            ];
    end
    xi = [xi0 zeta] - (df.' \ f.').';
    
    xi0 = xi(1:end-1);
    zeta = xi(end);
    
    if (niter > 100)
        break
    end
end % of NR iterations

xi = [xi0 zeta];

end % of function be_invmap
