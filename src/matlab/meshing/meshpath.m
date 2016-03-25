function m = meshpath(path, Le)
%MESHPATH Create line mesh from a path of lines and Bezier Curves
% MESH = MESHPATH(PATH, L) creates a NiHu mesh MESH from the eps path data
% PATH. The path is meshed by applying an approximate uniform element
% length L.
%
% See also: READ_EPSPATH

%   Copyright 2008-2010 P. Fiala
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last updated: 2015.03.07.

%
nP = size(path,1);
for iP = 1 : nP
    switch path(iP,1)
        case 1
            r0 = bezier(reshape(path(iP,[2 3 2 3 4 5 4 5]).',2,4).', Le);
        case 2
            r0 = bezier(reshape(path(iP,2:9).',2,4).', Le);
    end
    n = size(r0,1);
    m0 = create_empty_mesh();
    m0.Nodes = [(1:n).', r0, zeros(n,1)];
    m0.Elements = [(1:n-1).', repmat([ShapeSet.LinearLine.Id,1,1], n-1, 1), (1:n-1).' (2:n).'];
    if iP == 1
        m = m0;
    else
        m = join_meshes(m, m0);
    end
end

m.Elements(:,3) = 1;
m.Elements(:,4) = 1;
end

%
function r = bezier(R, Le)
%BEZIER Coordinates of Bezier curve
%  R = BEZIER(P, L) returns the coordinates of a cubic Bezier curve
%  defined by 4 control points. The control points are given in the rows of
%  the 4x2 or 4x3 matrix P. The curve is discretized so that the
%  approximate element length is L.
%
% Last updated: 21.11.2009
t = linspace(0, 1, 1001).'; % paramteter space
r = [(1-t).^3, 3*(1-t).^2.*t, 3*(1-t).*t.^2, t.^3] * R;  % curve coord
dr = diff(r,1,1);
ds = [0; cumsum(sqrt(dot(dr, dr, 2)))]; % element lengths
[ds, n] = unique(ds);
L = ds(end);        % length of the Bezier curve
N = ceil(L / Le);   % number of segments
l = (0:N)*L/N;      % divide into equal long segments
r = interp1(ds, r(n,:), l, 'linear', 'extrap');
end
