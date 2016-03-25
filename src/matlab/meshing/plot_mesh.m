function plot_mesh(mesh, varargin)
%PLOT_MESH Plot colored NiHu mesh
% PLOT_MESH(MESH) plots a NiHu mesh. 
% PLOT_MESH(MESH, 'node', DATA) plots a NiHu mesh colored by DATA
%    containing the color indices of the MESH's nodes.
% PLOT_MESH(MESH, 'elem', DATA) plots a NiHu mesh colored by DATA
%    containing the color indices of the MESH's elements.
% PLOT_MESH(MESH, DATA) plots a colored NiHU mesh. If DATA has the same
%    number of elements as MESH.Nodes, then the coloring is performed
%    nodewise. If DATA has the same number of elements as MESH.Elements,
%    then the coloring is performed elementwise. 
% PLOT_MESH(MESH, 'node', NodeID, DATA) plots a colored NiHU mesh based
%    on nodal color values. NodeID and DATA are Mx1 vectors containing
%    the colored nodes' IDs and color indices.
% PLOT_MESH(MESH, 'elem', NodeID, DATA) plots a colored NiHU mesh based
%    on element color values. ElemID and DATA are Mx1 vectors containing
%    the colored elements' IDs and color indices.
%
% Example:
%   sphere = create_sphere(1, 10);
%   x = sphere.Nodes(:,2);
%   nodind = mesh_select(sphere, 'x > 0', 'ind');
%   figure;
%   plot_mesh(sphere, 'node', sphere.Nodes(nodind,1), x(nodind));
%
%   x = sphere.Nodes(:,2);
%   elem = drop_IDs(sphere);
%   x = mean(x(elem(:,5:8)), 2);
%   [nodind, elind] = mesh_select(sphere, 'x < 0', 'ind');
%   figure;
%   plot_mesh(sphere, 'elem', sphere.Elements(elind,1), x(elind));

%   Copyright 2008-2010 P. Fiala
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last modified: 30.11.2009.

%% Argument check and default parameters
narginchk(1, 4);
nargoutchk(0, 0);

switch nargin
    case 1  % PLOT_MESH(MESH);
        cData = zeros(size(mesh.Nodes,1),1);
        cDOF = mesh.Nodes(:,1);
        nodewise = true;
    case 2  % PLOT_MESH(MESH, DATA);
        cData = varargin{1};
        if size(cData,1) == size(mesh.Nodes,1)
            cDOF = mesh.Nodes(:,1);
            nodewise = true;
        elseif size(cData,1) == size(mesh.Elements,1)
            cDOF = mesh.Elements(:,1);
            nodewise = false;
        else
            error('nihu:plot_mesh:argError', ...
                'Plot mode can not be determined from length(DATA).');
        end
    case 3
        switch lower(varargin{1})
            case 'node' % PLOT_MESH(MESH, 'node', DATA)
                nodewise = true;
                cDOF = mesh.Nodes(:,1);
            case 'elem' % PLOT_MESH(MESH, 'elem', DATA)
                nodewise = false;
                cDOF = mesh.Elements(:,1);
            otherwise
                error('nihu:plot_mesh:argError', ...
                    'Invalid type specifier ''%s''', varargin{1});
        end
        cData = varargin{2};
    case 4
        switch lower(varargin{1})
            case 'node'  % PLOT_MESH(MESH, 'node', DOF, DATA)
                nodewise = true;
            case 'elem'  % PLOT_MESH(MESH, 'elem', DOF, DATA)
                nodewise = false;
            otherwise
                error('nihu:plot_mesh:argError', ...
                    'Invalid type specifier ''%s''', varargin{1});
        end
        cDOF = varargin{2};
        cData = varargin{3};
end

if size(cData,1) ~= numel(cDOF)
    error('nihu:plot_mesh:argError', 'Invalid size of color data');
end

%% Convert ID-s to element and node indices, and drop IDs
if nodewise
    ID = mesh.Nodes(:,1);
else
    ID = mesh.Elements(:,1);
end
% check whether all IDs are contained in the mesh
d = ismember(cDOF, ID);
if ~all(d)
    error('nihu:plot_mesh:invalidID',...
        ['Invalid IDs (do not refer to mesh IDs):\n'...
        sprintf('%d\n', cDOF(~d))]);
end

dind(ID) = 1:length(ID);    % create ID-to-index permutation vector
c = zeros(size(ID,1), size(cData,2));        % null data
c(dind(cDOF(d)),:) = cData(d,:);

x = mesh.Nodes(:,2);
y = mesh.Nodes(:,3);
z = mesh.Nodes(:,4);
Elements = drop_IDs(mesh);


lsetid = Elements(:,2);

%% Plot LINE elements
lin = lsetid == ShapeSet.LinearLine.Id | lsetid == ShapeSet.InfiniteLine.Id;
if any(lin)
    elem = Elements(lin,5:6);
    line(x(elem.'), y(elem.'), z(elem.'), 'Color', 'black');
end

%% Plot TRIA elements
tria = (lsetid == ShapeSet.LinearTria.Id | lsetid == 223);
if any(tria)
    elem = Elements(tria,5:7);
    if nodewise
        patch('Faces', elem,...
            'Vertices', [x y z], ...
            'FaceVertexCData', c,...
            'FaceColor', 'interp');
    else
        patch('Faces', elem,...
            'Vertices', [x y z], ...
            'FaceVertexCData', c(tria),...
            'FaceColor', 'flat');
    end
end

%% Plot QUAD elements
quad = find(lsetid == ShapeSet.LinearQuad.Id |...
    lsetid == ShapeSet.InfiniteLinearQuad.Id | lsetid == 224);
if any(quad)
    elem = Elements(quad,5:8);
    if nodewise
        patch('Faces', elem,...
            'Vertices', [x y z], ...
            'FaceVertexCData', c,...
            'FaceColor', 'interp');
    else
        patch('Faces', elem,...
            'Vertices', [x y z], ...
            'FaceVertexCData', c(quad),...
            'FaceColor','flat');
    end
end

%% Plot 3D elements
% 3D elements are plotted by extracting their boundary and plotting the
% bounday's 2D surfrace elements with plot_mesh.
[uniqueIds, ~, idx] = unique(mesh.Elements(:,2));
uniqueLsets = ShapeSet.fromId(uniqueIds);
uniqueDims = zeros(size(uniqueLsets,1),1);
for e = 1 : length(uniqueLsets)
    uniqueDims(e) = uniqueLsets(e).Domain.Space.Dimension;
end
dims = uniqueDims(idx);
i3D = dims == 3;    % 3D selection
if any(i3D)
    mesh.Elements = mesh.Elements(i3D,:);   % keep 3D elements
    [bou, elemind] = get_boundary(mesh);      % compute boundary
    if nodewise
        plot_mesh(bou, 'node', cDOF, cData);
    else
        plot_mesh(bou, 'elem', bou.Elements(:,1), c(elemind));
    end
end

%% Axis settings
xlabel('x');
ylabel('y');
zlabel('z');
% axis equal tight;
grid on;
if max(abs(mesh.Nodes(:,4))) > 0
    view(3);
end

end
