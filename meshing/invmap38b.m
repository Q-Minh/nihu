function xi = invmap38b(xvert, x0)

x0 = x0(:);

maxiter = 10;
converged = false;
tol = 1e-5;

xivec = [
    0 0 0
    -1 -1 -1
    1 -1 -1
    1 1 -1
    -1 1 -1
    -1 -1 1
    1 -1 1
    1 1 1
    -1 1 1
    ];

for p = 1 : 9
    xi = xivec(p,:);
    iIter = 0;
    while true
        iIter = iIter+1;
        [N, dN] = shape(xi);
        %         [N, dN, ddN] = shape(xi);
        x = xvert' * N;
        d = x - x0;
        f = d' * d;
        if f < tol^2
            converged = true;
            break;
        end
        dx = xvert' * dN;
        df = (2 * d' * dx)';
        ddf = dx'* dx;
        %         for k = 1 : 8
        %             s = ddN(:,:,k);
        %             for i = 1 : 3
        %                 ddf = ddf + d(i) * xvert(k,i) * s;
        %             end
        %         end
        ddf = 2*ddf;
        
        xi = (xi' - ddf \ df)';
        if iIter > maxiter
            break;
        end
    end
    if converged
        break;
    end
end
% if ~converged
%     disp('nihu');
% end

end

function [N, dN, ddN] = shape(xivec)

xi = xivec(1);
eta = xivec(2);
zeta = xivec(3);

N = [
    (1-xi)*(1-eta)*(1-zeta)
    (1+xi)*(1-eta)*(1-zeta)
    (1+xi)*(1+eta)*(1-zeta)
    (1-xi)*(1+eta)*(1-zeta)
    (1-xi)*(1-eta)*(1+zeta)
    (1+xi)*(1-eta)*(1+zeta)
    (1+xi)*(1+eta)*(1+zeta)
    (1-xi)*(1+eta)*(1+zeta)
    ] / 8;

dN = [
    -(1-eta)*(1-zeta), -(1-xi)*(1-zeta), -(1-xi)*(1-eta)
    (1-eta)*(1-zeta), -(1+xi)*(1-zeta), -(1+xi)*(1-eta)
    (1+eta)*(1-zeta), (1+xi)*(1-zeta), -(1+xi)*(1+eta)
    -(1+eta)*(1-zeta), (1-xi)*(1-zeta), -(1-xi)*(1+eta)
    -(1-eta)*(1+zeta), -(1-xi)*(1+zeta), (1-xi)*(1-eta)
    (1-eta)*(1+zeta), -(1+xi)*(1+zeta), (1+xi)*(1-eta)
    (1+eta)*(1+zeta), (1+xi)*(1+zeta), (1+xi)*(1+eta)
    -(1+eta)*(1+zeta), (1-xi)*(1+zeta), (1-xi)*(1+eta)
    ] / 8;

if nargout > 2
    ddN(1,:,:) = [
        0, (1-zeta), (1-eta)
        0, -(1-zeta), -(1-eta)
        0, (1-zeta), -(1+eta)
        0, -(1-zeta), (1+eta)
        0, (1+zeta), -(1-eta)
        0, -(1+zeta), (1-eta)
        0, (1+zeta), (1+eta)
        0, -(1+zeta), -(1+eta)
        ]' / 8;
    ddN(2,:,:) = [
        (1-zeta), 0, (1-xi)
        -(1-zeta), 0, (1+xi)
        (1-zeta), 0, -(1+xi)
        -(1-zeta), 0, -(1-xi)
        (1+zeta), 0, -(1-xi)
        -(1+zeta), 0, -(1+xi)
        (1+zeta), 0, (1+xi)
        -(1+zeta), 0, (1-xi)
        ]' / 8;
    ddN(3,:,:) = [
        (1-eta), (1-xi), 0
        -(1-eta), (1+xi), 0
        -(1+eta), -(1+xi), 0
        (1+eta), -(1-xi), 0
        -(1-eta), -(1-xi), 0
        (1-eta), -(1+xi), 0
        (1+eta), (1+xi), 0
        -(1+eta), (1-xi), 0
        ]' / 8;
end
end