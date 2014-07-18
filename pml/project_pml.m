function pml_mesh = project_pml(surface,L,N,method,varargin)
% PROJECT_PML Project PML elements around a given surface
%   PML_MESH = PROJECT_PML(SURFACE,L,N,METHOD,VARARGIN) returns a mesh
%       consisting of PML elements projected from a surface mesh SURFACE
%       using the projection method selected by the parameter METHOD.
%       L is the thickness and N is the resolution of the layer.
%       VARARGIN should consists of the parameters
%       for the given method. METHOD must be selected from:
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
error(nargchk(5,6,nargin));
error(nargoutchk(0,1,nargout));

% Property IDs of PML elements
pID = max(surface.Properties(:,1))+1;
pID = pID:pID-1+N;

surface = drop_unused_nodes(surface);
Elements = drop_IDs(surface);
coords = surface.Nodes(:,2:4);
nNodes = size(surface.Nodes,1);
nElem = size(Elements,1);

r0 = varargin{1};
pml_lin = linspace(0,L,N+1);
pml_lin = reshape(repmat(pml_lin(2:end),nNodes,1),N*nNodes,1);
map_factor = ones(N*nNodes,3).*repmat(pml_lin,1,3);
ccoords = repmat(coords,N,1);

switch method
    case 'point'
        d = ccoords-repmat(r0,N*nNodes,1);
        map_coords = ccoords+map_factor.*d./repmat(sqrt(sum(d.^2,2)),1,3);
    case 'line'
        v = varargin{2}/norm(varargin{2});
        r1 = repmat(r0,N*nNodes,1)+repmat(dot(ccoords-repmat(r0,N*nNodes,1),repmat(v,N*nNodes,1),2),1,3).*repmat(v,N*nNodes,1);
        d = (ccoords-r1);
        map_coords = ccoords+map_factor.*d./repmat(sqrt(sum(d.^2,2)),1,3);
    case 'section'
        l = norm(varargin{2}-varargin{1});
        v = (varargin{2}-varargin{1})./l;
        dvec = dot(ccoords-repmat(r0,N*nNodes,1),repmat(v,N*nNodes,1),2);
        r1 = repmat(r0,N*nNodes,1)+repmat(max([zeros(N*nNodes,1),min([dvec,l*ones(N*nNodes,1)],[],2)],[],2),1,3).*repmat(v,N*nNodes,1);
        d = (ccoords-r1);
        map_coords = ccoords+map_factor.*d./repmat(sqrt(sum(d.^2,2)),1,3);
    case 'plane'
        n = varargin{2}/norm(varargin{2});
        r1 = repmat(n,N*nNodes,1).*repmat(dot(ccoords-repmat(r0,N*nNodes,1),repmat(n,N*nNodes,1),2),1,3);
        d = (ccoords-r1);
        map_coords = ccoords+map_factor.*d./repmat(sqrt(sum(d.^2,2)),1,3);
    case 'nodes'
        if(norm(size(r0)-size(coords))~=0)
            error('project_pml:ivalidArg','r0 must be the same size as nodes matrix for this method');
        end
        d = (ccoords-repmat(r0,N,1));
        map_coords = ccoords+map_factor.*d./repmat(sqrt(sum(d.^2,2)),1,3);
end

pml_mesh.Nodes = [(1:(N+1)*nNodes).',[coords;map_coords]];
pml_mesh.Elements = zeros(N*nElem,4+max(mod(Elements(:,2),10)));
o = ones(N,1);

% Creating PML elements from boundary elements
for cc = 1:nElem
    bn = mod(Elements(cc,2),10);    % number of element nodes
    nodes = Elements(cc,5:4+bn);
    
    if bn == 3      % TRIA element
        pat_nodes = repmat(nodes,N,2);
        pat = [(0:N-1).' (0:N-1).' (0:N-1).' (1:N).' (1:N).' (1:N).'];
    else            % Line element
        pat = [(0:N-1).' (1:N).' (1:N).' (0:N-1).'];
        if bn == 4  % QUAD element
            nodes = [nodes(1:2) nodes([4 3])];
            pat = repmat(pat,1,2);
        end
        pat_nodes = reshape(repmat(nodes,2*N,1),N,2*bn);
        
    end
    
    pml_mesh.Elements((cc-1)*N+1:cc*N,1:4+2*bn) = [(cc:cc+N-1).',(210+Elements(cc,2)+bn)*o,...
                                                   Elements(cc,3)*o,pID(:),pat_nodes+nNodes*pat];
end

pml_mesh.Materials = surface.Materials;
pml_mesh.Properties = surface.Properties;
pml_mesh.Properties(end+1:end+N,1:6) = [pID(:) 3*o o L*o N*o (1:N).'];
end