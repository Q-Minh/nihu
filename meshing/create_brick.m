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

%% Parameter check
switch nargin
    case 1 % Cx (= Cy = Cz) mode
        Cx = sort(varargin{1});
        Cy = Cx;
        Cz = Cx;
        N = length(Cx) * [1 1 1] - 1;
    case 2 % R N mode
        % Process corners
        R = varargin{1};
        switch size(R,1)
            case 1
                if isscalar(R)
                    R = [R R R];
                end
                R = [
                    0    0    0
                    R(1) 0    0
                    R(1) R(2) 0
                    0    R(2) 0
                    0    0    R(3)
                    R(1) 0    R(3)
                    R(1) R(2) R(3)
                    0    R(2) R(3)
                    ];
            case 2
                R = [
                    R(1,1) R(1,2) R(1,3)
                    R(2,1) R(1,2) R(1,3)
                    R(2,1) R(2,2) R(1,3)
                    R(1,1) R(2,2) R(1,3)
                    R(1,1) R(1,2) R(2,3)
                    R(2,1) R(1,2) R(2,3)
                    R(2,1) R(2,2) R(2,3)
                    R(1,1) R(2,2) R(2,3)
                    ];
            case 8
            otherwise
                error('NiHu:create_brick:argFormat',...
                    'Unsupported format of input arguments.');
        end
        N = varargin{2};
        if isscalar(N)
            N = [N N N];
        end
    case 3 % Cx Cy Cz mode
        Cx = sort(varargin{1});
        Cy = sort(varargin{2});
        Cz = sort(varargin{3});
        N = [length(Cx) length(Cy) length(Cz)] - 1;
    otherwise
        error('NiHu:create_brick:argNumber',...
            'Unsupported number of arguments: %d.', nargin);
end

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
if exist('R', 'var')
    phi = shapefun(model.Nodes(:,2:4), 38);
    model.Nodes(:,2:4) = phi * R;
elseif exist('Cx', 'var')
    [z y x] = meshgrid(Cz,Cy,Cx);            % create coordinates
    model.Nodes(:,2:4) = [x(:) y(:) z(:)];   % replace coordinates
end

end

