function [H, G] = bemHG_bm(model, k, alpha, pairs)
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
if nargin < 4
    pairs = [];
end

%% Gaussian quadrature divisions
%         singular near mid far
%gauss3 = assemble_gauss_struct(3, [9 7 5 2]);
gauss3 = assemble_gauss_struct(3, [1 13 7 2]);
gauss4 = assemble_gauss_struct(4, [9 7 5 2]);
[nodes, elements] = extract_bem_mesh(model);
dist = nodes(elements(:,3)+1,:) - nodes(elements(:,2)+1,:);
dist = [2 5] * max(sqrt(dot(dist,dist,2)));

% make sure alpha and k are complex
alpha = complex(real(alpha), imag(alpha));
k = complex(real(k), imag(k));

%[H, G] = bemHG_lin_bm(nodes, elements, gauss3, gauss4, dist, k, alpha);


if ~isempty(pairs)
    %% sparse matrices (FMBEM mode)
    [Hp, Gp] = bemHG_const_bm_sp(nodes, elements, gauss3, gauss4,...
        dist, k, alpha, pairs-1);
    n = size(elements,1);
    H = sparse(pairs(:,1), pairs(:,2), Hp, n, n);
    G = sparse(pairs(:,1), pairs(:,2), Gp, n, n);
else
    [H, G] = bemHG_const_bm(...
        nodes, elements, gauss3, gauss4, dist, k, alpha);
end

end
