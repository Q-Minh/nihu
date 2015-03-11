function [model, origind] = quad2tria(model)
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

%% Argument check
narginchk(1, 1);

%% 
q = model.Elements(:,2) == 24; % search for quad elements
quad = model.Elements(q,:);
elem = drop_IDs(model);
% Calculate length of diagonals
d1 = model.Nodes(elem(q,5),2:4) - model.Nodes(elem(q,7),2:4);
d1 = dot(d1, d1, 2);
d2 = model.Nodes(elem(q,6),2:4) - model.Nodes(elem(q,8),2:4);
d2 = dot(d2, d2, 2);

dd = d2 > d1;

tria = [                        % convert to TRIA
    quad(dd,[1:4, 5 6 7]);
    quad(dd,[1:4, 5 7 8]);
    quad(~dd,[1:4, 5 6 8]);
    quad(~dd,[1:4, 6 7 8]);
    ];
tria(:,2) = 23;
origind = [find(~q); repmat(find(q),2,1)];
model.Elements = model.Elements(~q,:);
model.Elements(end+(1:size(tria,1)),1:7) = tria;
model.Elements(:,1) = 1:size(model.Elements,1); % renumber elements
end