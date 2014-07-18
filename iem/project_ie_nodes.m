function [map_coords,ie_coords] = project_ie_nodes(coords,P,method,varargin)

%Check arguments
error(nargchk(4,5,nargin));
error(nargoutchk(0,2,nargout));

%Check P value
%TODO: Check P value

nNodes = size(coords,1);
z = linspace(-1+2/P,1-2/P,P-1);
if P>1
    factors = repmat(reshape(repmat(2./(1-z)-1,nNodes,1),(P-1)*nNodes,1),1,3);
    rep_coords = repmat(coords,P-1,1);
else
    ie_coords = [];
end

r0 = varargin{1};

map_factor = ones(nNodes,3);


switch method
    case 'point'
        if P>1
            ie_coords = rep_coords+factors.*(rep_coords-repmat(r0,nNodes*(P-1),1));
        end
        map_coords = coords+map_factor.*(coords-repmat(r0,nNodes,1));
    case 'line'
        v = varargin{2}/norm(varargin{2});
        r1 = repmat(r0,nNodes,1)+repmat(dot(coords-repmat(r0,nNodes,1),repmat(v,nNodes,1),2),1,3).*repmat(v,nNodes,1);
        if P>1
            ie_coords = rep_coords+factors.*(rep_coords-repmat(r1,P-1,1));
        end
        map_coords = coords+map_factor.*(coords-r1);
    case 'section'
        l = norm(varargin{2}-varargin{1}); % length of the section
        v = (varargin{2}-varargin{1})./l;  % direction vector of the section
        % distance from base point one
        dvec = dot(coords-repmat(r0,nNodes,1),repmat(v,nNodes,1),2);
        r1 = repmat(r0,nNodes,1) + repmat(max([zeros(nNodes,1),min([dvec,l*ones(nNodes,1)],[],2)],[],2),1,3).*repmat(v,nNodes,1);
        if P>1
            ie_coords = rep_coords+factors.*(rep_coords-repmat(r1,P-1,1));
        end
        map_coords = coords+map_factor.*(coords-r1);
    case 'plane'
        n = varargin{2}/norm(varargin{2});
        r1 = repmat(n,nNodes,1).*repmat(dot(coords-repmat(r0,nNodes,1),repmat(n,nNodes,1),2),1,3);
        if P>1
            ie_coords = rep_coords+factors.*(rep_coords-repmat(r1,P-1,1));
        end
        map_coords = coords+map_factor.*(coords-r1);
        
    case 'nodes'
        if(norm(size(r0)-size(coords))~=0)
            error('project_ie_nodes:ivalidArg','r0 must be the same size as nodes matrix for this method');
        end
        if P>1
            ie_coords = rep_coords+factors.*(rep_coords-repmat(r0,(P-1),1));
        end
        map_coords = coords+map_factor.*(coords-r0);
end