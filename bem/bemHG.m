function [H, G] = bemHG(model, k, type, points, pairs)
%BEMHG Build matrices of an acoustic bem model
%  [H, G] = bemHG(model, k, type, points, pairs) computes the bem system
%  matrices H and G.
%  For the case of surface matrices, the output matrices are defined with
%  Hp = Gp'
%  For the case of field point matrices, the output matrices are defined as
%  p_field = H p_surf - G p'_surf
% Parameters:
%   model  : AcouFEM model structure. The bem model is generated using the
%            QUAD and TRIA elements of the original model
%   k      : scalar wave number (real)
%   type   : 'lin'   for linear shape functions (one DOF per corner node)
%            'const' for constant shape functions (one DOF per element)
%   points : Nx3 xyz field point coordinates. If empty (default),
%            then the surface acoustic matrices are computed.
%   pairs  : Mx2 ij indices of the nonzero elements of the sparse output
%            matrices. If empty (default), the full bem matrices are computed.

% Peter Fiala
% Last modified: 17.11.2009.

%% parameter check
if nargin < 5
    pairs = [];
end
if nargin < 4
    points = [];
end

%% Gaussian quadrature divisions
%         singular near mid far
gauss3 = assemble_gauss_struct(3, [9 7 5 2]);
gauss4 = assemble_gauss_struct(4, [9 7 5 2]);
[nodes, elements] = extract_bem_mesh(model);
dist = nodes(elements(:,3)+1,:) - nodes(elements(:,2)+1,:);
dist = [2 5] * max(sqrt(dot(dist,dist,2)));

%% System matrices
if ~isempty(pairs)
    %% sparse matrices (FMBEM mode)
    if ~isempty(points) % field point sparse system matrices
        Nn = size(points,1);
        switch lower(type)
            case 'const'
                [Hp, Gp] = bemHG_const_sp(nodes, elements, gauss3, gauss4,...
                    dist, k, points, pairs-1);
                Ne = size(elements,1);
            case 'lin'
                error('The linear sparse option of bemHG is not yet implemented.');
        end
    else % surface sparse system matrices
        switch lower(type)
            case 'const'
                [Hp, Gp] = bemHG_const_sp(nodes, elements, gauss3, gauss4,...
                    dist, k, pairs-1);
                Nn = size(elements,1);
                Ne = Nn;
            case 'lin'
                error('The linear sparse option of bemHG is not yet implemented');
        end
    end
    H = sparse(pairs(:,1), pairs(:,2), Hp, Nn, Ne);
    G = sparse(pairs(:,1), pairs(:,2), Gp, Nn, Ne);
else
    %% full matrices (BEM mode)
    if ~isempty(points) % full field point matrices
        switch lower(type)
            case 'const'
                [H, G] = bemHG_const(nodes, elements, gauss3, gauss4,...
                    dist, k, points);
            case 'lin'
                [H, G] = bemHG_lin(nodes, elements, gauss3, gauss4,...
                    dist, k, points);
        end
    else % full surface matrices
        switch lower(type)
            case 'const'
                [H, G] = bemHG_const(nodes, elements, gauss3, gauss4, dist, k);
            case 'lin'
                [H, G] = bemHG_lin(nodes, elements, gauss3, gauss4, dist, k);
        end
    end
end
end
