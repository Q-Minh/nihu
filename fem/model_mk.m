function varargout = model_mk(model, opt)
%MODEL_MK Compute mass and stiffness matrices of the FE model
%   [M K C DOF] = MODEL_MK(MODEL) Computes the mass matrix M, stiffness
%   matrix K and damping matrix C of a FE/IE model MODEL. The output vector
%   DOF relates the elements of the output matrices to the degrees of
%   freedom of the model.
%   [M K DOF] = MODEL_MK(MODEL)
%   [I J M K C DOF] = MODEL_MK(MODEL, 'ind')
%   [I J M K DOF] = MODE_MK(MODEL, 'ind')
% See also: MODEL_MK_CYL, MODEL_A, MODEL_A_CYL

% Peter Fiala and Peter Rucz
% 2009 - 2011

% Argument check
narginchk(1,2);
nargoutchk(1,6);

% Process arguments
if nargin < 2
    % Choose matrix assembly as default
    opt = 'mat';
end
    

% Preprocessing
% Drop non-acoustical elements
AcMat = model.Materials(model.Materials(:,2) == 1,1);
model.Elements = model.Elements(ismember(model.Elements(:,3), AcMat),:);
% drop User IDs
Elements = drop_IDs(model);
% PML workaround
Elements = Elements(Elements(:,2) < 200,:);

%% Preallocating space for the acoustic matrices
% The elements of the sparse matrices are collected in arrays. The matrix
% values of each element are stored in non overlapping parts of the arrays.
% do not calculate infinite elements
nNode = sum(Elements(Elements(:,2) < 100,5:end)~=0, 2);   % number of nodes per element
inNode = sum(Elements(Elements(:,2) > 100,5:end)~=0,2)-mod(Elements(Elements(:,2) > 100,2),10);
nVal = nNode.^2;                        % number of matrix entries per elem
inVal = inNode.^2;
Start = [0; cumsum(nVal)];              % starting indices of elements
istart = [0; cumsum(inVal)];
N = Start(end); % total number of matrix entries (non overlapped)
iN = istart(end);
I  = zeros(N,1);  % preallocating row indices
iI = zeros(iN,1);
J  = zeros(N,1);  % preallocating column indices
iJ = zeros(iN,1);
M  = zeros(N,1);  % preallocating mass values
K  = zeros(N,1);  % preallocating stiffness values

iM = zeros(iN,1);
iK = zeros(iN,1);
iC = zeros(iN,1);

%% Matrix assembly 
nElem = size(Elements,1);               % number of elements
iElem = 0;                              % 

for type = [12 23 24 34 36 38]
    typeind = Elements(:,2) == type;
    if ~any(typeind)
        continue;
    end
    dim = floor(type/10);               % dimension number
    nnodes = mod(type,10);              % number of nodes of an element
    elements = Elements(typeind,:);     % the elements of this type
    nelem = size(elements,1);           % number of elements of this type
    start = Start(typeind);             % start indices of elements
    [w, N, dN] = elem_quad(type);       % Gaussian quadarture for this type
    for ielem = 1 : nelem
        elem = elements(ielem,:);
        nodes = elem(4+(1:nnodes));             % elem nodes
        vert = model.Nodes(nodes,1+(1:dim)).';  % corner vertices
        % Assemble element matrix
        [m, k] = elem_mk(type, vert, w, N, dN, model.Materials(elem(3),3:4));
        i = repmat(nodes.', 1, nnodes);         % row indices
        j = i.';                                % column indices
        ind = start(ielem)+(1:nnodes^2);        % index positions
        I(ind) = i(:);
        J(ind) = j(:);
        M(ind) = m(:);
        K(ind) = k(:);
        iElem = iElem+1;
        progbar(1, nElem, iElem);
    end
end

iecount = 0;
for type = [111 122 133 134]
    typeind = Elements(:,2) == type;
    if ~any(typeind)
        continue;
    end
    
    telements = Elements(typeind,:);          % the elements of this type
    % select infinite elements by properties
    props = unique(Elements(typeind,4));      % the unique properties
    dim = floor((type-100)/10);               % dimension number
    nb = mod(type,10);                        % number of base points
    for prop = props
        propind = telements(:,4) == prop;
        elements = telements(propind,:);
        
        P = model.Properties(prop,3);         % polynomial degree of elements
        ptype = model.Properties(prop,4);     % polynomial type
        param = model.Properties(prop,5:6);   % polynom parameters
        
        % Load the corresponding gauss file
        G = ielem_loadgauss(type,ptype,P,param);
        
        % indices of mapping nodes
        minds = 5:(4+2*nb);
        % indices of pressure nodes
        pinds = [5:(4+nb), (5+2*nb):(4+(P+1)*nb)];
                
        nelem = size(elements,1);           % number of elements of this type
        for ielem = 1: nelem
            iecount = iecount +1;           % increase counter
            elem = elements(ielem,:);       % the current element
            %Mapping matrix
            map = model.Nodes(elem(minds),2:2+dim-1);
                        
            % Assemble the element matrices
            [m, k, c] = ielem_mkc(dim,nb,map,...
                model.Materials(elem(3),3:4),P,G,...
                model.Properties(elem(4),5:6));
            
            i = repmat(elem(pinds).', 1, inNode(iecount));  % row indices
            j = i.';                                        % col indices
            
            ind = istart(iecount)+1:istart(iecount+1);
            % Fill index arrays
            iI(ind) = i(:);
            iJ(ind) = j(:);
            iM(ind) = m(:);
            iK(ind) = k(:);
            iC(ind) = c(:);
            iElem = iElem+1;
            progbar(1, nElem, iElem);
        end
    end
end
%% Check for infinite elements
C = [zeros(size(I)); iC];
% Concatenate index arrays
I = [I; iI];
J = [J; iJ];
M = [M; iM];
K = [K; iK];
%% Generate the sparse matrices from the arrays
nNodes = size(model.Nodes,1);               % Number of nodes
sel = ismember(1:nNodes,[I;iI]);            % Logical array of selected nodes
vNodes = sum(sel);                          % Number of value nodes
index = zeros(vNodes,1);                    % Indices
index(sel) = 1:vNodes;
% The DOF vector contains the node IDs corresponding to matrix rows
DOF = model.Nodes(sel,1);
% Renumber indices according to valid nodes
I = index(I);
J = index(J);

% Switch between output options
switch opt
    case {'mat','matrix'}
        M = sparse(I, J, M, vNodes, vNodes);
        K = sparse(I, J, K, vNodes, vNodes);
        varargout{1} = M;
        varargout{2} = K;
        switch nargout
            case 3 % [M K DOF]
                varargout{3} = DOF;
            case 4 % [M K C DOF]
                C = sparse(I, J, C, vNodes, vNodes);
                varargout{3} = C;
                varargout{4} = DOF;
        end
    case {'ind','index'}
        varargout{1} = I;
        varargout{2} = J;
        varargout{3} = M;
        varargout{4} = K;
        switch nargout
            case 5 % [I J M K DOF]
                varargout{5} = DOF;
            case 6 % [I J M K C DOF]
                varargout{5} = C;
                varargout{6} = DOF;
        end
                
end

end % of function model_mk


function [w, N, dNxi] = elem_quad(type)

%% Compute Gaussian quadrature weights and shape function samples
nnod = mod(type, 10);
switch type
    case 12
        [xi, w] = gaussquad1(2, nnod);
    case {23, 24}
        [xi, w] = gaussquad2(2, nnod);
    case 34
        [xi, w] = gaussquad3(2, nnod);
    case {36, 38}
        [xi, w] = gaussquad3(3, nnod);
end
[N, dNxi] = shapefun(xi, type);
dNxi = reshape(shiftdim(dNxi,2),numel(dNxi)/nnod,nnod);
end

function [M, K] = elem_mk(type, xnode, w, N, dNxi, material)

dim = floor(type/10);
ngauss = length(w);

%% Coordinate transform
J = dNxi * xnode.';      % Jacobian of the coordiate transform
dNx = zeros(size(dNxi)); % Gradient in xyz space
j = zeros(ngauss,1);
for n = 1 : ngauss
    ind = dim*(n-1)+(1:dim);
    Jac = J(ind,:);
    dNx(ind,:) = (Jac) \ dNxi(ind,:);
    j(n) = abs(det(Jac));
end

%% Numerical integration with Gaussian quadrature
q = w.*j;
M = N.'*diag(q)*N * material(1);
q = reshape(repmat(q.',dim, 1),dim*ngauss,1);
K = dNx.' * diag(q) * dNx * (material(1)*material(2)^2);
end