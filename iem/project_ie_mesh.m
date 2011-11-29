function ie_mesh = project_ie_mesh(surface,P,method,varargin)
% PROJECT_IE_MESH Project infinite elements from a given surface
%   IE_MESH = PROJECT_IE_MESH(SURFACE,P,METHOD,VARARGIN) returns a mesh
%       consisting of infinite elements projected from a surface mesh
%       SURFACE of polynomial order P using the projection method selected
%       by the parameter METHOD. VARARGIN should consists of the parameters
%       for the given method. The projection methods are the same as for 
%       project_ie_nodes function. METHOD must be selected from:
%   'point'    - Projection from a point given by its coordinates.
%                Requires one additional parameter, R0, a 1x3 vector with
%                the (x, y, z) coordinates of the projection point.
%   'line'     - Projection from a line
%                Requires two additional parameters, R0, a point on the
%                projection line and V, the direction vector of this line.
%   'section'  - Projection from a line section
%                Requires two additional parameters, R0 and R1 the two
%                endpoints of the section.
%   'plane'    - Projection from a plane
%                Requires two additional parameters, R0, a point on the
%                projection plane and N, the normal vector of this plane.
%   'nodes'    - Projection from independent nodes
%                Requires one additional parameter, C, an Nx3 matrix with
%                the coordinates of the independent projection nodes. N
%                must equal the number of nodes of the surface.

% Check arguments
error(nargchk(4,5,nargin));
error(nargoutchk(0,2,nargout));

surface = drop_unused_nodes(surface);
Elements = drop_IDs(surface);
coords = surface.Nodes(:,2:4);
nNodes = size(surface.Nodes,1);
[map_coords,ie_coords] = project_ie_nodes(coords,P,method,varargin{:});
ie_mesh.Nodes = [(1:(P+1)*nNodes).',[coords;map_coords;ie_coords]];
nElem = size(Elements,1);
ie_mesh.Elements = zeros(nElem,4+(P+1)*max(mod(Elements(:,2),10)));
for c1 = 1:nElem
    bn = mod(Elements(c1,2),10); %number of boundary element nodes
    ie_mesh.Elements(c1,1:4+(P+1)*bn) = [c1,Elements(c1,2)+110,Elements(c1,3:4),...
        repmat(Elements(c1,5:4+bn),1,P+1)+nNodes*reshape(repmat(0:P,bn,1),1,(P+1)*bn)];
end
ie_mesh.Materials = surface.Materials;
ie_mesh.Properties = surface.Properties;