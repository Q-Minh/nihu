function [model origind] = quad2tria(model)
%QUAD2TRIA  Replace Quad elements with TRIA elements
%   [MESH ORIGID] = QUAD2TRIA(MESH) replaces each QUAD element in MESH by
%   two TRIA elements. The element IDs are not conserved. The vector ORIGID
%   contains the original element ID of each element in the output mesh.
%
% Example:
%  s = create_sphere_boundary(1,10);
%  x = s.Nodes(s.Elements(:,5),2);  % approx. x coordinate of an element
%  [s2, i] = quad2tria(s);
%  plot_mesh(s2, 'elem', x(i))

%   Copyright 2008-2010 P. Fiala
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last modified: 02.12.2009

%% Argument check
error(nargchk(1, 1, nargin, 'struct'));

%% 
q = model.Elements(:,2) == 24; % search for quad elements
quad = model.Elements(q,:);
tria = [                        % convert to TRIA
    quad(:,[1:4, 5 6 7]);
    quad(:,[1:4, 5 7 8]);
    ];
tria(:,2) = 23;
origind = [find(~q); repmat(find(q),2,1)];
model.Elements = model.Elements(~q,:);
model.Elements(end+(1:size(tria,1)),1:7) = tria;
model.Elements(:,1) = 1:size(model.Elements,1); % renumber elements
end
