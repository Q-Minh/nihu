function bb = boundingbox(mesh)
%BOUNDINGBOX Compute bounding box of a NiHu mesh
%  BB = BOUNDINGBOX(MESH) returns a 2x3 matrix containing the corner
%  coordinates of the bounding box of the NiHu mesh MESH.

%   Copyright 2008-2010 P. Fiala
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last modified: 02.12.2009
%% Parameter check
error(nargchk(1, 1, nargin,'struct'));

%%
bb = [
    min(mesh.Nodes(:,2:4), [], 1);
    max(mesh.Nodes(:,2:4), [], 1)
    ];

end
