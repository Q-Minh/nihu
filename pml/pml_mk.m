function varargout = pml_mk(model, omega, opt)
%PML_MK Compute mass and stiffness matrices of the PML model
%   [M_pml K_pml DOF] = PML_MK(MODEL, OMEGA) Computes the mass matrix M_pml, stiffness
%   matrix K_pml of a PML model MODEL. OMEGA is the angular frequency.
%   The output vector DOF relates the elements of the output matrices to the degrees of
%   freedom of the model.
%   [M K DOF] = PML_MK(MODEL)
%   [I J M K DOF] = PML_MK(MODEL, 'ind') 

%% Argument check
error(nargchk(2,3,nargin,'struct'));
error(nargoutchk(1,5,nargout, 'struct'));

% Process arguments
if nargin < 3
    % Choose matrix assembly as default
    opt = 'mat';
end 

%% Preprocessing
% Drop non-acoustical elements
AcMat = model.Materials(model.Materials(:,2) == 1,1);
model.Elements = model.Elements(ismember(model.Elements(:,3), AcMat),:);
% drop User IDs
Elements = drop_IDs(model);
Elements = Elements(Elements(:,2) > 200,:);

%% Preallocating space for the acoustic matrices
% The elements of the sparse matrices are collected in arrays. The matrix
% values of each element are stored in non overlapping parts of the arrays.
nNode = sum(Elements(Elements(:,2) > 200,5:end)~=0, 2);   % number of nodes per element
nVal = nNode.^2;                        % number of matrix entries per elem
Start = [0; cumsum(nVal)];              % starting indices of elements
N = Start(end); % total number of matrix entries (non overlapped)
I  = zeros(N,1);  % preallocating row indices
J  = zeros(N,1);  % preallocating column indices
M  = zeros(N,1);  % preallocating mass values
K  = zeros(N,1);  % preallocating stiffness values

%% Matrix assembly 
nElem = size(Elements,1);               % number of elements
iElem = 0;                              % 

for type = [212 223 224 234 236 238]
    typeind = Elements(:,2) == type;
    if ~any(typeind)
        continue;
    end
    
    dim = floor((type-200)/10);                        % dimension number
    nnodes = mod(type-200,10);                         % number of nodes of an element
    elements = Elements(typeind,:);                    % the elements of this type
    nelem = size(elements,1);                          % number of elements of this type
    start = Start(typeind);                            % start indices of elements
    propid = model.Properties(elements(1,4),3);
    hsigma = model.PMLData{propid,2};                  % Absorbing function handle
    
    % Compute Gaussian quadrature weights and shape function samples
    switch type
        case 212
            [xi, w] = gaussquad1(3, nnodes);
        case {223, 224}
            [xi, w] = gaussquad2(3, nnodes);
        case 234
            [xi, w] = gaussquad3(3, nnodes);
        case {236, 238}
            [xi, w] = gaussquad3(4, nnodes);
    end
    [N, dN] = shapefun(xi, type-200);
    dN = reshape(shiftdim(dN,2),numel(dN)/nnodes,nnodes);
    
    if type == 236
        xi = 2*xi(:,3)-1;       % transform zeta to (-1,1)
    else
        xi = xi(:,1);
    end
    
    for ielem = 1:nelem
        elem = elements(ielem,:);
        nodes = elem(4+(1:nnodes));             % elem nodes
        vert = model.Nodes(nodes,1+(1:dim));    % corner vertices
        
        % Assemble element matrix
        if model.Properties(elem(4),6) == 0;
            propid = model.Properties(elem(4),3);
            hsigma = model.PMLData(propid,2:4);
            psigma = model.PMLData(propid,5:7);
            [m_pml, k_pml] = pmlelem_mk_glob(type, vert, w, N, dN, omega, hsigma, psigma, model.Properties(elem(4),4), model.Materials(elem(3),3:4));
        else
            [m_pml, k_pml] = pmlelem_mk_loc(dim, vert, w, N, dN, xi, omega, hsigma, model.Properties(elem(4),4:6), model.Materials(elem(3),3:4));
        end
        i = repmat(nodes.', 1, nnodes);         % row indices
        j = i.';                                % column indices
        ind = start(ielem)+(1:nnodes^2);        % index positions
        I(ind) = i(:);
        J(ind) = j(:);
        M(ind) = m_pml(:);
        K(ind) = k_pml(:);
        iElem = iElem+1;
        progbar(1, nElem, iElem);
    end
end
    
%% Generate the sparse matrices from the arrays
nNodes = size(model.Nodes,1);
m = ismember(1:nNodes,I);
vNodes = sum(m);
index = zeros(vNodes,1);
index(m) = 1:vNodes;
% The DOF vector contains the node IDs corresponding to matrix rows 
DOF = model.Nodes(m,1);
I = index(I);
J = index(J);

% Switch between output options
switch opt
    case {'mat','matrix'}
        M = sparse(I, J, M, vNodes, vNodes);
        K = sparse(I, J, K, vNodes, vNodes);
        varargout{1} = M;
        varargout{2} = K;
        varargout{3} = DOF;
    case {'ind','index'}
        varargout{1} = I;
        varargout{2} = J;
        varargout{3} = M;
        varargout{4} = K;
        varargout{5} = DOF;
end

end