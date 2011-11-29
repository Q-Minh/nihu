function [N, dN] = shapefun(x, type)
%SHAPEFUN Shape functions over standard elements
% [N, dN] = SHAPEFUN(X, TYPE) returns linear shape function samples over
% standard quad and tria elements. The input vector X is a nx2 vector, each
% row containins the coordinates [xi eta] of a point in the standard
% element. The input parameter TYPE is 23 for TRIA and 24 for QUAD
% elements. The output N contains the n samples of the shape function
% N(xi, eta), the output parameter dN contains the derivatives
% [N'_xi N'_eta] for each location.

%   Copyright 2008-2010 P. Fiala
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last modified 02.12.2009
%% Argument check
error(nargchk(2, 3, nargin, 'struct'));

%%
switch type
    case 12
        %% Line element
        xi = x(:,1);
        N = [1-xi, 1+xi]/2;
        dN = repmat([-1, +1]/2, size(xi,1),1);
    case 23
        %% TRIA element
        xi = x(:,1);
        eta = x(:,2);
        N = [1-xi-eta, xi, eta];
        dN(:,:,1) = repmat([-1, +1,  0], size(xi,1),1);
        dN(:,:,2) = repmat([-1,  0, +1], size(xi,1),1);
    case 24
        %% QUAD element
        xi = x(:,1);
        eta = x(:,2);
        N = [(1-xi).*(1-eta), (1+xi).*(1-eta),...
            (1+xi).*(1+eta), (1-xi).*(1+eta)] / 4;
        dN(:,:,1) = [-(1-eta), (1-eta), (1+eta), -(1+eta)] / 4;
        dN(:,:,2) = [-(1-xi), -(1+xi), (1+xi), (1-xi)] / 4;
    case 34
        %% TETRA element
        xi = x(:,1);
        eta = x(:,2);
        zeta = x(:,3);
        N = [1-xi-eta-zeta, xi, eta, zeta];
        dN(:,:,1) = repmat([-1, 1, 0, 0], size(xi,1),1);
        dN(:,:,2) = repmat([-1, 0, 1, 0], size(xi,1),1);
        dN(:,:,3) = repmat([-1, 0, 0, 1], size(xi,1),1);
    case 36
        %% PENTA element
        xi = x(:,1);
        eta = x(:,2);
        zeta = x(:,3);
        n = size(xi,1);
        N = [
            (1-xi-eta).*(1-zeta), xi.*(1-zeta), eta.*(1-zeta), ...
            (1-xi-eta).*(1+zeta), xi.*(1+zeta), eta.*(1+zeta)
            ]/2;
        dN(:,:,1) = [-(1-zeta), (1-zeta), zeros(n,1),...
            -(1+zeta), (1+zeta), zeros(n,1)]/2;
        dN(:,:,2) = [-(1-zeta), zeros(n,1), (1-zeta),...
            -(1+zeta), zeros(n,1), (1+zeta)]/2;
        dN(:,:,3) = [-(1-xi-eta), -xi, -eta,...
            (1-xi-eta), xi, eta]/2;
    case 38
        %% HEXA element
        xi = x(:,1);
        eta = x(:,2);
        zeta = x(:,3);
        N = [
            (1-xi).*(1-eta).*(1-zeta), (1+xi).*(1-eta).*(1-zeta),...
            (1+xi).*(1+eta).*(1-zeta), (1-xi).*(1+eta).*(1-zeta),...
            (1-xi).*(1-eta).*(1+zeta), (1+xi).*(1-eta).*(1+zeta),...
            (1+xi).*(1+eta).*(1+zeta), (1-xi).*(1+eta).*(1+zeta)
            ]/8;
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
    otherwise
        error('NiHu:shapefun:argValue',...
            'Elem type %d not supported', type);
end
end
