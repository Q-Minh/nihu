function model = create_slab(varargin)
%CREATE_SLAB  Create slab mesh (NiHu / meshing)
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

%   Copyright 2008-2012 P. Fiala, P. Rucz
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last modifed: 2012.12.12.

% extract arguments
args = create_slab_args(varargin{:});
N = args.N;

% Check the processed variables
if size(N,1) ~= 1 || size(N,2) ~= 2
    error('NiHu:create_slab:argFormat',...
        'Unsupported format of input arguments.');
end

if any(N < 1)
    error('NiHu:create_slab:argValue',...
        'Number of segments (%d %d) must be positive on each side', N(1), N(2));
end

% Create base model
model = create_slab_base(N);
% Apply transformation
if isfield(args, 'R')
    phi = shapefun(model.Nodes(:,2:3), 24);
    model.Nodes(:,2:4) = phi * args.R;
elseif isfield(args, 'Cx')
    [y x] = meshgrid(args.Cy,args.Cx);            % create coordinates
    model.Nodes(:,2:3) = [x(:) y(:)];   % replace coordinates
end

end
