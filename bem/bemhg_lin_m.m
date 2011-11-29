function [H, G] = bemhg_lin_m(nodes, elements, gauss3, gauss4, dist, k)
%BEMHG_LIN_M Generate acoustic BEM system matrices

nNode = size(nodes,1);

H = zeros(nNode, nNode);
G = zeros(nNode, nNode);

nElement = size(elements,1);
for iNode = 1 : nNode
    q = nodes(iNode,:);
    for iElement = 1 : nElement
        nVert = elements(iElement,1);
        elem = elements(iElement,1+(1:nVert))+1;
        vert = nodes(elem,:);

        [sing, st] = ismember(iNode, elem);
        if sing  % singular integral
            ind = mod(st+(0:nVert-1)-1, nVert)+1;
            elem = elem(ind);
            switch nVert
                case 3
                    [P, dP] = int_reg(q, vert(ind,:), gauss3(1), k);
                case 4
                    [P, dP] = int_reg(q, vert(ind,:), gauss4(1), k);
            end
        else % regular integral
            % select gaussian quadrature
            d = norm(mean(vert,1) - q);
            if d < dist(1)
                iG = 1;
            elseif d < dist(2)
                iG = 2;
            else
                iG = 3;
            end
            switch nVert
                case 3
                    [P, dP] = int_reg(q, vert, gauss3(iG), k);
                case 4
                    [P, dP] = int_reg(q, vert, gauss4(iG), k);
            end
        end
        % perform integral
        H(iNode,elem) = H(iNode,elem) + dP;
        G(iNode,elem) = G(iNode,elem) + P;
    end
end
H = H - 2*pi*eye(size(H));


function [P, dP] = int_reg(q, vert, g, k)

nVert = size(vert, 1);
gcoord = g.N * vert;
gnorm = cross(g.Nxi * vert, g.Neta * vert, 2);
jac = sqrt(dot(gnorm, gnorm, 2));
gnorm = gnorm ./ repmat(jac, 1,3);
[p, dp] = incident('point', q, gcoord, gnorm, k);
% p = ones(size(p));
% dp = ones(size(dp));
dP = (g.w .* jac).' * (repmat(dp, 1, nVert) .* g.N);
P = (g.w .* jac).' * (repmat(p, 1, nVert) .* g.N);
