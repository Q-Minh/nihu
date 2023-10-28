function plot_node_numbers(model, shift)
%PLOT_NODE_NUMBERS Plot node numbers of a NiHu mesh
%   PLOT_NODE_NUMBERS(MESH, SHIFT) plots the node numbers of a NiHu mesh in
%   the current axes. The node numbers are shifted from the node
%   coordinates by an amount defined in the 3D vector SHIFT. The default of
%   SHIFT is [0 0 0].

%   Copyright 2008-2010 P. Fiala
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last modified 02.12.2009
%% Argument check
error(nargchk(1, 2, nargin, 'struct'));
if nargin == 1
    shift = [0 0 0];
end
shift = shift(:)';
if length(shift) < 3
    shift(3) = 0;
end

%%
coord = bsxfun(@plus, model.Nodes(:,2:4), shift);
text(coord(:,1), coord(:,2), coord(:,3), num2str(model.Nodes(:,1)));
end
