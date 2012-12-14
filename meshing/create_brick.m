function model = create_brick(varargin)
%CREATE_BRICK Create a brick volume mesh (NiHu / meshing)
%   BRICK = CREATE_BRICK(L, N) creates a brick model with side lengths
%   given in the 3D row vector L and division given in the 3D row vector N.
%   The brick is located at the origin of the coordinate system, its faces
%   are aligned along the coordinate axes. N contains the number of brick
%   elements along the three directions.
%
%   For a scalar input L, it is assumed that the side lengths are the
%   same along the three dimensions.
%   For a scalar N, the same is assumed for the number of elements.
%
%   If L has two rows then the brick is located between L(1,:) and L(2,:).
%
%   BRICK = CREATE_BRICK(C, N) where C is a 8x3 matrix creates a brick with
%   given corner nodes defined in the rows of the 8x3 matrix C.
%
%   BRICK = CREATE_BRICK(Cx, Cy, Cz) where Ci are column vectors creates a
%   brick whose nodes are elements of the Descartes product Cx x Cy x Cz.
%   If only Cx is defined, it is asumed that Cy = Cz = Cx.
%
%   Example:
%      L = [5 4 3];
%      Le = .3;
%      N = ceil(L / Le);
%      brick1 = create_brick(L, N);
%
%      C = [
%           0   0  0
%           1   0  0
%           1.5 1  0
%           -.1 .9 .1
%           0   0  1
%           1   0  1.2
%           1.5 1  1
%           -.1 .7 1
%         ];
%      brick2 = create_brick(C, [10, 10, 10]);
%
%      Cx = (.1:.05:1).^2;
%      Cy = (.1:.05:1).^3;
%      Cz = (.1:.05:1).^4;
%      brick3 = create_brick(Cx, Cy, Cz);
%
% See also: CREATE_BRICK_BOUNDARY, CREATE_LINE, CREATE_SLAB, CREATE_CIRCLE,
% CREATE_CIRCLE_QUADRANT, CREATE_SPHERE, CREATE_SPHERE_BOUNDARY,
% CREATE_CATSEYE

%   Copyright 2008-2012 P. Fiala, P. Rucz
%   Budapest University of Technology and Economics

% Last modifed: 2012.12.14.

% extract arguments
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

