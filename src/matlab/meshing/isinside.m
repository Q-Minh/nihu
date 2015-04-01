function inside = isinside(mesh, points)
%ISINSIDE  Check whether points are inside a mesh volume
%   INSIDE = ISINSIDE(MESH, POINTS) returns a Nx1 logical array INSIDE
%   whose elements indicate whether a point given in rows of the Nx3
%   matrix POINTS is inside the NiHu mesh MESH or not.
%
% See also: CHIEFGEN

%   Copyright 2008-2015 P. Fiala
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last modified: 2015.03.31. PF: new ShapeSet system introduced

narginchk(2,2);

% Initialization
nPoints = size(points,1);   % number of points to check
inside = false(nPoints, 1); % initializing indicators

if any(mesh.Elements(:,2) == ShapeSet.LinearLine.Id)
    dim = 2;                % mesh is 2D
else
    dim = 3;                % mesh is 3D
end

switch dim
    case 2
        % Check for polygon
        elem = drop_IDs(mesh);
        x1 = mesh.Nodes(elem(:,5),2);
        y1 = mesh.Nodes(elem(:,5),3);
        x2 = mesh.Nodes(elem(:,6),2);
        y2 = mesh.Nodes(elem(:,6),3);
        dir = (x2-x1) ./ (y2-y1);
        for iPoints = 1 : nPoints     %check separately for each point
            y = points(iPoints,2);
            g = find((y1 > y) ~= (y2 > y)); % elements around y
            inside(iPoints) = mod(...
                sum((points(iPoints,1)-x1(g)) < (y-y1(g)).*dir(g)),...
                2)==1;
        end
    case 3
        % Check for polyhedron
        elem = drop_IDs(quad2tria(mesh));
        xe = mesh.Nodes(:,2);
        xe = xe(elem(:,5:7));   % x coordinates of elements
        xe = [min(xe,[],2) max(xe,[],2)];
        ye = mesh.Nodes(:,3);
        ye = ye(elem(:,5:7));   % y coordinates of elements
        ye = [min(ye,[],2) max(ye,[],2)];
        for iPoints = 1 : nPoints     %check separately for each point
            x = points(iPoints,1:3);    % point coordiantes
            g = find(x(1) < xe(:,2) & x(1) > xe(:,1));
            g = g(x(2) < ye(g,2) & x(2) > ye(g,1));
            c = false;  % indicates parity of the intersections
            for iG = 1 : length(g)
                e = elem(g(iG),5:7);
                xx = mesh.Nodes(e,2:4).';
                m = [xx(:,2)-xx(:,1), xx(:,3)-xx(:,1), [0;0;1]];
                a = m \ (x.'-xx(:,1));
                c = xor(c, all(a > 0) && a(1)+a(2) <= 1);
            end
            inside(iPoints) = c;
        end
end

end
