function model = create_brick(varargin)
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

%   Copyright 2008-2012 P. Fiala, P. Rucz
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last modifed: 2012.12.13.

args = create_brick_args(varargin{:});
N = args.N;

% Check the processed variables
if size(N,1) ~= 1 || size(N,2) ~= 3
    error('NiHu:create_brick:argFormat',...
        'Unsupported format of input arguments.');
end

if any(N < 1)
    error('NiHu:create_brick:argValue',...
        'Number of segments (%d %d %d) must be positive on each side', N(1), N(2));
end

% Create base model
model = create_brick_base(N);
% Apply transformation
if isfield(args, 'R')
    phi = shapefun(model.Nodes(:,2:4), 38);
    model.Nodes(:,2:4) = phi * args.R;
elseif isfield(args, 'Cx')
    [z y x] = meshgrid(args.Cz,args.Cy,args.Cx);    % create coordinates
    model.Nodes(:,2:4) = [x(:) y(:) z(:)];          % replace coordinates
end

end

