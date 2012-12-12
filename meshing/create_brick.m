function brick = create_brick(varargin)
%CREATE_BRICK  Create a brick volume mesh
%   BRICK = CREATE_BRICK(L, N) creates a brick model with side lengths
%   given in the 3D vector L and division given in the 3D vector N.
%   The brick is located at the origin of the coordinate system, its faces
%   are aligned along the coordinate axes. N contains the number of brick
%   elements along the three directions. For a scalar input L, it is
%   assumed that the side lengths are the same along the three dimensions.
%   For a scalar N, the same is assumed for the number of elements.
%
%   BRICK = CREATE_BRICK(C, N) where C is a 8x3 matrix creates a brick with
%   given corner nodes defined in the rows of the matrix C. N is a scalar
%   or a 3D vector containing the number of elements along the three
%   directions.
%
%   BRICK = CREATE_BRICK(Cx, Cy, Cz) where Ci are column vectors creates a
%   brick whose nodes are elements of the Descartes product Cx x Cy x Cz.
%   If only Cx is defined, it is asumed that Cy = Cz = Cx.
%
% See also: CREATE_BRICK_BOUNDARY, CREATE_LINE, CREATE_SLAB, CREATE_CIRCLE,
% CREATE_CIRCLE_QUADRANT, CREATE_SPHERE, CREATE_SPHERE_BOUNDARY,
% CREATE_CATSEYE

%   Copyright 2008-2010 P. Fiala
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last modifed: 02.12.2009

%% Parameter check
if isscalar(L)
    L = repmat(L,1,3);
end
if nargin == 1
    N = L;
    Cz = L;
end
if isscalar(N)
    N = repmat(N,1,3);
end

%% Creating the brick model
if numel(N) > 3 % Node vectors are given
    Cx = L;
    Cy = N;
    N = [length(Cx) length(Cy) length(Cz)]-1;
    brick = create_brick(N, N);
    brick.Nodes(:,2) = Cx(brick.Nodes(:,2)+1);
    brick.Nodes(:,3) = Cy(brick.Nodes(:,3)+1);
    brick.Nodes(:,4) = Cz(brick.Nodes(:,4)+1);
elseif size(L, 1) == 8  % corner nodes are given
    % uniform centered slab
    brick = translate_mesh(create_brick(2, N), -[1 1 1]);
    % node transformation
    phi = shapefun(brick.Nodes(:,2:4), 38);
    brick.Nodes(:,2:4) = phi * L;
else    % side lengths are given
    brick = extrude_mesh(create_slab(L(1:2), N(1:2)), [0 0 L(3)/N(3)], N(3));
end

end

