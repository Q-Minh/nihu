function field = field_bbx(bb, Le, N)
%FIELD_BBX Cross-shaped field point mesh based on bounding box data
%  FIELD = FIELD_BBX(BB, Le, N) returns a Cross shaped field point mesh
%  stretching the bounding box BB and meshed with element length Le. The
%  optional input variable N determines the number of cross slabs in the
%  three dimensions.
%
% Example:
%  field = field_bbx([0 0 0; 5 4 3], .1, [3 1 1]);
%  plot_mesh(field);

%   Copyright 2008-2010 P. Fiala
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last modified 08.12.2009.

%% Argument check
error(nargchk(2,3,nargin));
if nargin < 3
    N = [1 1 1];
end
N = max(N, 0);

%% Dimensions of bounding box
D = diff(bb,1,1);   % side length
C = mean(bb,1);     % center

%% Create slabs in three directions
field.Nodes = zeros(0,4);
field.Elements = zeros(0,4);
field.Materials = [];
field.Properties = [];
for d = 1 : 3
    if N(3) >= 1
        field0 = create_slab(D(1:2), ceil(D(1:2)/Le));
        field0 = repeat_mesh(field0, [0, 0, D(3)/(N(3))], N(3)-1);
        field0 = translate_mesh(field0, -mean(field0.Nodes(:,2:4,1)));
        field = join_meshes(field, field0);
    end
    D = circshift(D, [0 1]);
    N = circshift(N, [0 1]);
    field = rotate_mesh(field, [0 0 0], [0 0 1], pi/2);
    field = rotate_mesh(field, [0 0 0], [0 1 0], pi/2);
end
field = translate_mesh(field, C);

end
