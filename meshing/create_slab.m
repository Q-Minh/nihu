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

switch nargin
    case 1 % One argument mode (Cx, Cx mode)
        Cx = sort(varargin{1});
        Cy = Cx;
        N = length(Cx) * [1 1] - 1;
    case 2 % Two argument mode
        % Cx, Cy mode
        if (size(varargin{1},1) > 2 && size(varargin{1},2) == 1) || size(varargin{2},1) > 2
            Cx = sort(varargin{1});
            Cy = sort(varargin{2});
            N = [length(Cx)-1 length(Cy)-1];
        % R, N mode
        else
            % Process corners
            R = varargin{1};
            switch size(R,1)
                case 1
                    if isscalar(R)
                        R = [R R];
                    end
                    R = [
                        0    0    0
                        R(1) 0    0
                        R(1) R(2) 0
                        0    R(2) 0
                        ];
                case 2
                    R = [
                        R(1,1) R(1,2) 0
                        R(2,1) R(1,2) 0
                        R(2,1) R(2,2) 0
                        R(1,1) R(2,2) 0
                        ];
                case 4
                    R = [R zeros(size(R,1), 3-size(R,2))];
                otherwise
                    error('NiHu:create_slab:argFormat',...
                        'Unsupported format of input arguments.');
            end
            N = varargin{2};
            if isscalar(N)
                 N = [N N];
            end
        end
    otherwise
        error('NiHu:create_slab:argNumber',...
            'Unsupported number of arguments: %d.', nargin);
end

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
if exist('R', 'var')
    phi = shapefun(model.Nodes(:,2:3), 24);
    model.Nodes(:,2:4) = phi * R;
elseif exist('Cx', 'var')
    [y x] = meshgrid(Cy,Cx);            % create coordinates
    model.Nodes(:,2:3) = [x(:) y(:)];   % replace coordinates
end

end
