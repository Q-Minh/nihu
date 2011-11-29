function slab = create_slab(L, N)
%CREATE_SLAB  Create slab mesh
%   SLAB = CREATE_SLAB(L, N) Creates a slab with side lengths given in the 
%   2dim vector L and division parameters given in the 2dim vector N. If L
%   ore N are scalars, the brick's side lengths or division paramters
%   are assumed to be uniform. The slab is located at the origin, its sides
%   are aligned along the coordinate axes.
%
%   SLAB = CREATE_SLAB(C, N) Creates a slab with corner coordinates given
%   in the 4x3 matrix C and with division parameters given in the scalar or
%   2dim vector N.
%
%   SLAB = CREATE_SLAB(Cx, Cy) Creates an NxM slab using the N+1-element
%   vector Cx and the M+1-element vector Cy. The nodes of the slab are
%   elements of the Descartes product of Cx and Cy.
%
% See also: CREATE_LINE, CREATE_CIRCLE, CREATE_CIRCLE_QUADRANT, 
% CREATE_BRICK, CREATE_BRICK_BOUNDARY, CREATE_SPHERE,
% CREATE_SPHERE_BOUNDARY, CREATE_CATSEYE

%   Copyright 2008-2010 P. Fiala
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last modifed: 02.12.2009

%% Parameter check
error(nargchk(1, 2, nargin, 'struct'));

if numel(L) == 1
    L = repmat(L, 1, 2);
end
if nargin == 1
    N = L;
end
if numel(N) == 1
    N = repmat(N, 1, 2);
end

%% Create slab
if length(N) > 2 % Cx and Cy vectors are given
    Cx = sort(L);
    Cy = sort(N);
    N = [length(Cx) length(Cy)]-1;
    slab = create_slab(N, N);
    slab.Nodes(:,2) = Cx(slab.Nodes(:,2)+1);
    slab.Nodes(:,3) = Cy(slab.Nodes(:,3)+1);
elseif size(L,1) == 4   % Corner nodes are given
    % uniform centered slab
    slab = translate_mesh(create_slab(2, N), -[1 1 0]);
    % transformation
    xi = slab.Nodes(:,2);
    eta = slab.Nodes(:,3);
    phi = [(1-xi).*(1-eta), (1+xi).*(1-eta),...
        (1+xi).*(1+eta), (1-xi).*(1+eta)]/4;
    slab.Nodes(:,1+(1:size(L,2))) = phi * L;
else    % Side lengths are given
    slab = extrude_mesh(create_line(L(1), N(1)), [0 L(2)/N(2) 0], N(2));
end

end
