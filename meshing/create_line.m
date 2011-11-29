function model = create_line(L, N)
%CREATE_LINE  Create a line mesh
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

% Last modifed: 07.09.2010

%% Parameter check
error(nargchk(1, 2, nargin, 'struct'));

if nargin == 1
    N = length(L)-1;
    model.Nodes = [(1:N+1).' sort(L(:)), repmat([0 0], N+1, 1)];
elseif numel(L) == 1    % Side length is given
    model.Nodes = [(1:N+1).' linspace(0,L,N+1).', repmat([0 0], N+1, 1)];
elseif size(L,1) == 2
    model.Nodes = [(1:N+1).', ...
        linspace(L(1,1), L(2,1), N+1).', ...
        linspace(L(1,2), L(2,2), N+1).', ...
        linspace(L(1,3), L(2,3), N+1).'];
end

%%
model.Elements = [(1:N).', repmat([12 1 1], N, 1), (1:N).', (2:N+1).'];
model.Materials = [1 1];
model.Properties = [1 1];

end
