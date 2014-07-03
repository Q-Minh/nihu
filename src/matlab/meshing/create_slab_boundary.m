function model = create_slab_boundary(varargin)
%CREATE_SLAB_BOUNDARY Create a slab boundary mesh (NiHu / meshing)
%   SLAB = CREATE_SLAB_BOUNDARY(L, N) creates a slab model with side lengths
%   given in the 3D vector L and division given in the 3D vector N.
%   The slab is located at the origin of the coordinate system, its faces
%   are aligned along the coordinate axes. N contains the number of slab
%   elements along the three directions. For a scalar input L, it is
%   assumed that the side lengths are the same along the three dimensions.
%   For a scalar N, the same is assumed for the number of elements.
%
%   SLAB = CREATE_SLAB_BOUNDARY(C, N) where C is a 8x3 matrix creates a slab with
%   given corner nodes defined in the rows of the matrix C. N is a scalar
%   or a 3D vector containing the number of elements along the three
%   directions.
%
%   SLAB = CREATE_SLAB_BOUNDARY(Cx, Cy, Cz) where Ci are column vectors creates a
%   slab whose nodes are elements of the Descartes product Cx x Cy x Cz.
%   If only Cx is defined, it is asumed that Cy = Cz = Cx.
%
% See also: CREATE_SLAB_BOUNDARY, CREATE_LINE, CREATE_SLAB, CREATE_CIRCLE,
% CREATE_CIRCLE_QUADRANT, CREATE_SPHERE, CREATE_SPHERE_BOUNDARY,
% CREATE_CATSEYE

%   Copyright 2008-2012 P. Fiala, P. Rucz
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last modifed: 2012.12.13.

% extract arguments
args = create_slab_args(varargin{:});
N = args.N;

% Check the processed variables
if size(N,1) ~= 1 || size(N,2) ~= 3
    error('NiHu:create_slab_boundary:argFormat',...
        'Unsupported format of input arguments.');
end

if any(N < 1)
    error('NiHu:create_slab_boundary:argValue',...
        'Number of segments (%d %d %d) must be positive on each side', N(1), N(2));
end

% Create base model
model = create_slab_boundary_base(N);
% Apply transformation
if isfield(args, 'R')
    phi = shapefun(model.Nodes(:,2:4), 24);
    model.Nodes(:,2:4) = phi * args.R;
elseif isfield(args, 'Cx')
    error('NiHu:create_slab_boundary:argFormat',...
        'Cx Cy parametrisation is not supperted (yet)');
end

end
