function points = chiefgen(mesh, nChief)
%CHIEFGEN Generate random CHIEF points to a bem mesh
%  POINTS = CHIEFGEN(MESH, N) generates N CHIEF points to a
%  boundary element MESH. The CHIEF points are located inside the meshed
%  domain, and can be used as overdetermination points when the problem
%  of fictitious eigenfrequencies is solved by means of the CHIEF method.
%
% Example with full mesh:
%   cat = create_catseye(1,10);
%   p = chiefgen(cat, 20);
%   plot3(p(:,1), p(:,2), p(:,3), '.');
%   plot_mesh(cat); alpha .5;
%
% See also: bemHG

% Last modified: 02.12.2009
% Peter Fiala

%% Argument check
error(nargchk(1,2,nargin,'struct'));
if nargin == 1
    nChief = 50;
end

%%
mesh = drop_unused_nodes(mesh);
lim = [                             % bounding box of the mesh
    min(mesh.Nodes(:,2:4), [], 1)
    max(mesh.Nodes(:,2:4), [], 1)
    ];
l = diff(lim, 1, 1);                % dimensions of the bounding box

points = zeros(0,3);    % initializing CHIEF points
while size(points,1) < nChief
    % generate random points inside the bounding box
    xi = repmat(lim(1,:), nChief, 1) + rand(nChief,3) * diag(l);
    xi = xi(isinside(mesh, xi), :);
    % add new points to the Chief set
    points = [
        points
        xi(1:min(size(xi,1), nChief-size(points,1)), :)
        ]; %#ok<AGROW>
end
end