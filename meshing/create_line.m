function model = create_line(varargin)
%CREATE_LINE  Create a line mesh (NiHu / meshing)
%   LINE = CREATE_LINE(L, N) creates a line model with length given
%   in L and division given in N. The line starts at the origin, and is
%   aligned along the x axis.
%
%   LINE = CREATE_LINE(R, N) creates a line model with end locations given
%   in the 2x3 matrix R and division number given in N.
%
%   LINE = CREATE_LINE(Cx) creates a line model aligned along the x axis,
%   with nodes located at the positions given in vector Cx.
%
% See also: CREATE_SLAB, CREATE_CIRCLE, CREATE_CIRCLE_QUADRANT,
% CREATE_BRICK, CREATE_BRICK_BOUNDARY, CREATE_SPHERE,
% CREATE_SPHERE_BOUNDARY, CREATE_CATSEYE

%   Copyright 2008-2010 P. Fiala, P. Rucz
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last modifed: 2012.12.12

% TODO: Cx mode not working yet

% Process input arguments
switch nargin
    case 1 % one argument mode
        Cx = varargin{1};
        if size(Cx,2) > 3
            error('NiHu:create_line:argFormat',...
                'Unsupported format of input arguments.');
        end
        N = size(Cx,1) - 1;
    case 2 % Two arguments mode
        if (isscalar(varargin{2})) % L and N are given
            s = size(varargin{1});
            R = [
                zeros(2-s(1), 3);
                varargin{1}, zeros(s(1), 3-s(2))
                ];
            N = varargin{2};
        else
            error('NiHu:create_line:argFormat',...
                'Unsupported format of input arguments.');
        end
    otherwise
        error('NiHu:create_line:argNumber',...
            'Unsupported number of arguments: %d.', nargin);
end
if N < 1
    error('NiHu:create_line:argValue',...
        'Number of line segments (%d) must be positive', N);
end

% Transformation
model = create_line_base(N);            % Create base
if exist('R', 'var')
    phi = shapefun(model.Nodes(:,2), 12);   % Obtain shape values
    model.Nodes(:,2:4) = phi * R;           % Finish transformation
elseif exist('Cx', 'var')
    model.Nodes(:,1+(1:size(Cx,2))) = Cx;   % replace coordinates
end

end
