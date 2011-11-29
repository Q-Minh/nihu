function [El, Xi] = fe_invmap(mesh, points, tol)
%FE_INVMAP Inverse mapping matrix based on FE mesh and field points
%  [EL, XI] = FE_INVMAP(MESH, POINTS, TOL) returns the vector EL and the
%  Nxdim matrix XI that contains the element numbers and the corresponding
%  (xi, eta, zeta) local coordinates, respectively. Parameter TOL is the
%  tolerance for the norm of the distance between the field points and the
%  corresponding mapped locations.

% Peter Fiala, Peter Rucz
% Last modified: 14.04.2010.

%% Argument check and default parameters
error(nargchk(2,3,nargin,'struct'));
error(nargoutchk(2,2,nargout,'struct'));
if (nargin < 3) 
    tol = 1e-3;
end

%% Initialization
nElem = size(mesh.Elements,1);      % number of elements
nPoint = size(points,1);            % number of field points
dim = floor(mod(mesh.Elements(1,2),100)/10); % dimensionality of the problem
El = zeros(nPoint,1);
Xi = zeros(nPoint,dim);
elements = drop_IDs(mesh);  % get rid of node IDs
% Mesh coordinates
x = mesh.Nodes(:,2).';        % x coordinates of each vertex
y = mesh.Nodes(:,3).';        % y coordinates of each vertex
z = mesh.Nodes(:,4).';        % y coordinates of each vertex

%% Bounding box of each element
bb = zeros(nElem,6);
for n = [2 3 4 6 8]
    %Find the finite elements of the appropriate type
    q = find((elements(:,2) < 100)  & (mod(elements(:,2),10) == n));
    if ~isempty(q)
        e = elements(q, 4+(1:n));
        bb(q,:) = [min(x(e),[],2), min(y(e),[],2), min(z(e),[],2), ...
            max(x(e),[],2), max(y(e),[],2), max(z(e),[],2)];
    end
    %Find the infinite elements of the appropriate type
    iq = find((elements(:,2) > 100)  & (mod(elements(:,2),10) == n));
    if ~isempty(iq)
        for c1 = 1:numel(iq)
            Map = mesh.Nodes(elements(iq(c1),4+(1:2*n)),2:4);
            bb(iq(c1),:) = ieboundbox(Map);
        end
    end
end



%% Interpolation vector for each field point
for iPoint = 1 : nPoint
    xf = points(iPoint,1:dim); % field point coord
    % search for close elements whose bounding box contains the field point
    close = 1:nElem;
    for i = 1 : dim
        close = close(xf(i) >= bb(close,i) & xf(i) <= bb(close,i+3));
    end
    % solve for mapped coordinate in each close element
    for iC = 1 : length(close)
        elem = elements(close(iC),:); % close elements
        type = elem(2);
        nnod = mod(type,10);
        switch type
            case {12 23 24 36 38}
                nodes = elem(4+(1:nnod));
                xnod = mesh.Nodes(nodes,1+(1:dim));
            case {111 122 133 134}
                nodes = elem(4+(1:2*nnod));
                xnod = mesh.Nodes(nodes,1+(1:dim));
        end
        xi = invmap(xnod, xf, type);
        switch type
            case {12, 24, 38, 122, 134}
                xi = xi(all(abs(xi) < 1 + tol, 2), :);
                if isempty(xi)
                    continue
                end
            case {23, 34}   % TRIA and TETRA condition
                if any(abs(xi) < -tol) || sum(xi,2) > 1+tol
                    continue;
                end
        end
        El(iPoint) = close(iC);
        Xi(iPoint,:) = xi;
        break;
    end
    progbar(1, nPoint, iPoint);
end
end