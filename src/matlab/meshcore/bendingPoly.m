function [N, vars] = bendingPoly(c)

x = sym('x', 'real');
y = sym('y', 'real');

switch size(c,1)
    case 3
        L1 = 1-x-y;
        L2 = x;
        L3 = y;
        psi = [L1.^3, L2.^3, L3.^3,...
            L1.^2.*L2, L2.^2.*L3, L3.^2*L1,...
            L1.^2.*L3, L2.^2*L1, L3.^2.*L2];
    case 4
        % cubic polynomial base of shape functions
        psi = [1 x y x^2 x*y y^2 x^3 x^2*y x*y^2 y^3 x^3*y x*y^3];
    otherwise
        error('invalid number of nodal coordinates');
end

w = psi;            % vertical displacement
betax = diff(w,y);  % rotation around x axis
betay = -diff(w,x); % rotation around y axis

% nodal positions
nNodes = size(c,1);
nDofPerNode = 3;

G = zeros(nNodes*nDofPerNode,nNodes*nDofPerNode);
for k = 1 : nNodes
    G((k-1)*nDofPerNode+1,:) = subs(w, {x, y}, {c(k,1), c(k,2)});
    G((k-1)*nDofPerNode+2,:) = subs(betax, {x, y}, {c(k,1), c(k,2)});
    G((k-1)*nDofPerNode+3,:) = subs(betay, {x, y}, {c(k,1), c(k,2)});
end
N = simplify(psi / G);
vars = [x y];

end
