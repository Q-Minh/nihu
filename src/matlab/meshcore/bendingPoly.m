function [N, vars] = bendingPoly(domain)

d = domain.Space.Dimension;
g = sym('xi', [1 d]);

switch domain
    case Domain.Line
        x = g(1);
        w = [1 x x.^2 x.^3];
    case Domain.Tria
        x = g(1);
        y = g(2);
        L1 = 1-x-y;
        L2 = x;
        L3 = y;
        w = [L1.^3, L2.^3, L3.^3,...
            L1.^2.*L2, L2.^2.*L3, L3.^2*L1,...
            L1.^2.*L3, L2.^2*L1, L3.^2.*L2];
    case Domain.Quad
        % cubic polynomial base of shape functions
        x = g(1);
        y = g(2);
        w = [1 x y x^2 x*y y^2 x^3 x^2*y x*y^2 y^3 x^3*y x*y^3];
    otherwise
        error('invalid number of nodal coordinates');
end

switch d
    case 1
        U = [ w; diff(w,x) ];
    case 2
        U = [ w; diff(w,y); -diff(w,x) ];
end

% nodal positions
c = domain.CornerNodes;
nNodes = size(c,1);
nDofPerNode = 1 + d;

G = zeros(nNodes*nDofPerNode,nNodes*nDofPerNode);
for k = 1 : nNodes
    for s = 1 : nDofPerNode
        G((k-1)*nDofPerNode+s,:) = subs(U(s,:), g, c(k,:));
    end
end
N = simplify(w / G);
vars = g;

end
