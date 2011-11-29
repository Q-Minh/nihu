function B = ibemB(model, k, type, points, pairs)
%BEMHG Build matrices of an acoustic bem model
%  B = ibemB(model, k, type, points, pairs) computes the bem system
%  matrix B.
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
% singular near mid far
gauss3 = assemble_gauss_struct(3, [9 7 5 2]);
gauss4 = assemble_gauss_struct(4, [9 7 5 2]);
[nodes, elements] = extract_bem_mesh(model);
dist = nodes(elements(:,3)+1,:) - nodes(elements(:,2)+1,:);
dist = [2 5] * max(sqrt(dot(dist,dist,2)));

%% System matrices
if ~isempty(pairs)
    error('The sparse option of ibemB is not yet implemented.');
else
    %% full matrices (BEM mode)
    if ~isempty(points) % full field point matrices
        switch lower(type)
            case 'const'
                [~, B] = bemHG_const(nodes, elements, gauss3, gauss4, dist, k, points);
            case 'lin'
                [~, B] = bemHG_lin(nodes, elements, gauss3, gauss4, dist, k, points);
        end
    else % full surface matrices
        switch lower(type)
            case 'const'
                B = ibemB_const(nodes, elements, gauss3, gauss4, dist, k);
            case 'lin'
                B = ibemB_lin(nodes, elements, gauss3, gauss4, dist, k);
        end
    end
end
end