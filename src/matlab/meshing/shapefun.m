function [N, dN, ddN] = shapefun(x, type)
%SHAPEFUN Shape functions over standard elements
% [N, dN] = SHAPEFUN(X, TYPE) returns linear shape function samples over
% standard elements. The input vector X is a nxd vector, each
% row containins the coordinates [xi eta] of a point in the standard
% element. The output N contains the n samples of the shape function
% N(xi, eta), the output parameter dN contains the derivatives
% [N'_xi N'_eta] for each location.

%   Copyright 2008-2012 P. Fiala and P. Rucz
%   Budapest University of Technology and Economics

% Last modified 2013.12.16.

switch type
    case {12 121}
        % Line element
        xi = x(:,1);
        N = [1-xi, 1+xi]/2;
        if nargout > 1
            dN = repmat([-1, +1]/2, size(xi,1),1);
        end
        if nargout > 2
            ddN = zeros(length(xi), 2);
        end
    case 122
        % Quadratic Line element
        xi = x(:,1);
        N = [(xi-1).*xi/2, 1-xi.^2, (1+xi).*xi/2];
        if nargout > 1
            dN = [xi-.5, -2*xi, xi+.5];
        end
        if nargout > 2
            ddN = repmat([1, -2, 1], size(xi,1),1);
        end
    case {23 231}
        % TRIA element
        xi = x(:,1);
        eta = x(:,2);
        N = [1-xi-eta, xi, eta];
        if nargout > 1
            dN(:,:,1) = repmat([-1, +1,  0], size(xi,1),1);
            dN(:,:,2) = repmat([-1,  0, +1], size(xi,1),1);
        end
        if nargout > 2
            ddN = zeros(length(xi), 3, 3);
        end
    case {21}
        % constant QUAD shape
        N = ones(size(x,1), 1);
        if nargout > 1
            dN = zeros(size(x,1), 1, 2);
        end
        if nargout > 2
            ddN = zeros(size(x,1), 1, 3);
        end
    case {24 241}
        % linear QUAD element
        xi = x(:,1);
        eta = x(:,2);
        N = [(1-xi).*(1-eta), (1+xi).*(1-eta),...
            (1+xi).*(1+eta), (1-xi).*(1+eta)] / 4;
        if nargout > 1
            dN(:,:,1) = [-(1-eta), (1-eta), (1+eta), -(1+eta)] / 4;
            dN(:,:,2) = [-(1-xi), -(1+xi), (1+xi), (1-xi)] / 4;
        end
        if nargout > 2
            ddN(:,:,1) = zeros(length(xi), 4);
            ddN(:,:,2) = repmat([1, -1, +1, -1] / 4, length(xi), 1);
            ddN(:,:,3) = zeros(length(xi), 4);
        end
    case 232
        % QUADRATIC TRIA element
        xi = x(:,1);
        eta = x(:,2);
        N = [ (eta + xi - 1).*(2.*eta + 2.*xi - 1), -4.*xi.*(eta + xi - 1), xi.*(2.*xi - 1), 4.*eta.*xi, eta.*(2.*eta - 1), -4.*eta.*(eta + xi - 1)];
        if nargout > 1
            dN(:,:,1) = [4.*eta+4.*xi - 3, 4-8.*xi-4.*eta,       4.*xi - 1, 4.*eta, zeros(size(xi)),           -4.*eta];
            dN(:,:,2) = [4.*eta+4.*xi - 3,         -4.*xi, zeros(size(xi)),  4.*xi, 4.*eta - 1, 4 - 4.*xi - 8.*eta];
        end
    case 242
%         % QUADRATIC QUAD element
%         xi = x(:,1);
%         eta = x(:,2);
%         N = [ (eta.*xi.*(eta - 1).*(xi - 1))/4, -(eta.*(xi.^2 - 1).*(eta - 1))/2, (eta.*xi.*(eta - 1).*(xi + 1))/4, -(xi.*(eta.^2 - 1).*(xi + 1))/2, (eta.*xi.*(eta + 1).*(xi + 1))/4, -(eta.*(xi.^2 - 1).*(eta + 1))/2, (eta.*xi.*(eta + 1).*(xi - 1))/4, -(xi.*(eta.^2 - 1).*(xi - 1))/2, (eta.^2 - 1).*(xi.^2 - 1)];
%         if nargout > 1
%             dN(:,:,1) = [ (eta.*(2.*xi - 1).*(eta - 1))/4, -eta.*xi.*(eta - 1), (eta.*(2.*xi + 1).*(eta - 1))/4, -((eta.^2 - 1).*(2.*xi + 1))/2, (eta.*(2.*xi + 1).*(eta + 1))/4, -eta.*xi.*(eta + 1), (eta.*(2.*xi - 1).*(eta + 1))/4, -((eta.^2 - 1).*(2.*xi - 1))/2, 2.*xi.*(eta.^2 - 1)];
%             dN(:,:,2) = [ (xi.*(2.*eta - 1).*(xi - 1))/4, -((2.*eta - 1).*(xi.^2 - 1))/2, (xi.*(2.*eta - 1).*(xi + 1))/4, -eta.*xi.*(xi + 1), (xi.*(2.*eta + 1).*(xi + 1))/4, -((2.*eta + 1).*(xi.^2 - 1))/2, (xi.*(2.*eta + 1).*(xi - 1))/4, -eta.*xi.*(xi - 1), 2.*eta.*(xi.^2 - 1)];
%         end
        % QUADRATIC QUAD element
        xi = x(:,1);
        eta = x(:,2);
        N = [ -((xi - 1).*(eta - 1).*(xi + eta + 1))/4, ((xi.^2 - 1).*(eta - 1))/2, ((xi + 1).*(eta - 1).*(eta - xi + 1))/4, -((eta.^2 - 1).*(xi + 1))/2, ((xi + 1).*(eta + 1).*(xi + eta - 1))/4, -((xi.^2 - 1).*(eta + 1))/2, ((xi - 1).*(eta + 1).*(xi - eta + 1))/4, ((eta.^2 - 1).*(xi - 1))/2];
        if nargout > 1
            dN(:,:,1) = [ -((2*xi + eta).*(eta - 1))/4, xi*(eta - 1), -((2*xi - eta).*(eta - 1))/4, 1/2 - eta^2/2, ((2*xi + eta).*(eta + 1))/4, -xi*(eta + 1), ((2*xi - eta).*(eta + 1))/4, eta.^2/2 - 1/2];
            dN(:,:,2) = [ -((xi + 2*eta).*(xi - 1))/4, xi^2/2 - 1/2, -((xi - 2*eta).*(xi + 1))/4, -eta*(xi + 1), ((xi + 2*eta).*(xi + 1))/4, 1/2 - xi.^2/2, ((xi - 2*eta).*(xi - 1))/4, eta.*(xi - 1)];
        end
    case 34
        % TETRA element
        xi = x(:,1);
        eta = x(:,2);
        zeta = x(:,3);
        N = [1-xi-eta-zeta, xi, eta, zeta];
        if nargout > 1
            dN(:,:,1) = repmat([-1, 1, 0, 0], size(xi,1),1);
            dN(:,:,2) = repmat([-1, 0, 1, 0], size(xi,1),1);
            dN(:,:,3) = repmat([-1, 0, 0, 1], size(xi,1),1);
        end
    case 36
        % PENTA element
        xi = x(:,1);
        eta = x(:,2);
        zeta = x(:,3);
        n = size(xi,1);
        N = [
            (1-xi-eta).*(1-zeta), xi.*(1-zeta), eta.*(1-zeta), ...
            (1-xi-eta).*(1+zeta), xi.*(1+zeta), eta.*(1+zeta)
            ]/2;
        if nargout > 1
            dN(:,:,1) = [-(1-zeta), (1-zeta), zeros(n,1),...
                -(1+zeta), (1+zeta), zeros(n,1)]/2;
            dN(:,:,2) = [-(1-zeta), zeros(n,1), (1-zeta),...
                -(1+zeta), zeros(n,1), (1+zeta)]/2;
            dN(:,:,3) = [-(1-xi-eta), -xi, -eta,...
                (1-xi-eta), xi, eta]/2;
        end
    case 38
        % HEXA element
        xi = x(:,1);
        eta = x(:,2);
        zeta = x(:,3);
        N = [
            (1-xi).*(1-eta).*(1-zeta), (1+xi).*(1-eta).*(1-zeta),...
            (1+xi).*(1+eta).*(1-zeta), (1-xi).*(1+eta).*(1-zeta),...
            (1-xi).*(1-eta).*(1+zeta), (1+xi).*(1-eta).*(1+zeta),...
            (1+xi).*(1+eta).*(1+zeta), (1-xi).*(1+eta).*(1+zeta)
            ]/8;
        if nargout > 1
            dN(:,:,1) = [
                -(1-eta).*(1-zeta), (1-eta).*(1-zeta),...
                (1+eta).*(1-zeta), -(1+eta).*(1-zeta),...
                -(1-eta).*(1+zeta), (1-eta).*(1+zeta),...
                (1+eta).*(1+zeta), -(1+eta).*(1+zeta)
                ]/8;
            dN(:,:,2) = [
                -(1-xi).*(1-zeta), -(1+xi).*(1-zeta),...
                (1+xi).*(1-zeta), (1-xi).*(1-zeta),...
                -(1-xi).*(1+zeta), -(1+xi).*(1+zeta),...
                (1+xi).*(1+zeta), (1-xi).*(1+zeta)
                ]/8;
            dN(:,:,3) = [
                -(1-xi).*(1-eta), -(1+xi).*(1-eta),...
                -(1+xi).*(1+eta), -(1-xi).*(1+eta),...
                (1-xi).*(1-eta), (1+xi).*(1-eta),...
                (1+xi).*(1+eta), (1-xi).*(1+eta)
                ]/8;
        end
    otherwise
        error('NiHu:shapefun:argValue',...
            'Elem type %d not supported', type);
end
end
