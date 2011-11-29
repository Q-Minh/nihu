function [M, K, C, DOF] = iemodel_mk2(model)
%MODEL_MK Compute mass and stiffness matrices of the FE model
%   [M K DOF] = MODEL_MK(MODEL) Computes the mass matrix M and stiffness
%   matrix K of a FE model MODEL. The output vector DOF relates the
%   elements of the output matrices to the degrees of freedom of the model.
%
% See also: elem_mk, elem_mk_cyl, elem_a, elem_a_cyl, model_mk_cyl,
% model_a, model_a_cyl

% Peter Fiala
% 2009

%% Preprocessing
% Drop non-acoustical elements
AcMat = model.Materials(model.Materials(:,2) == 1,1);
model.Elements = model.Elements(ismember(model.Elements(:,3), AcMat),:);
% drop User IDs
Elements = drop_IDs(model);

%% Preallocating space for the acoustic matrices
% The elements of the sparse matrices are collected in arrays. The matrix
% values of each element are stored in non overlapping parts of the arrays.
% do not calculate infinite elements
nNode = sum(Elements(Elements(:,2) < 100,5:end)~=0, 2);   % number of nodes per element
inNode = sum(Elements(Elements(:,2) > 100,5:end)~=0,2)-mod(Elements(Elements(:,2) > 100,2),10);
nVal = nNode.^2;                        % number of matrix entries per elem
inVal = inNode.^2;
start = [0; cumsum(nVal)];              % starting indices of elements
istart = [0; cumsum(inVal)];
N = start(end); % total number of matrix entries (non overlapped)
iN = istart(end);
I = zeros(N,1); % preallocating row indices
iI = zeros(iN,1);
J = zeros(N,1); % preallocating column indices
iJ = zeros(iN,1);
M = zeros(N,1); % preallocating mass values
K = zeros(N,1); % preallocating stiffness values

iM = zeros(iN,1);
iK = zeros(iN,1);
iC = zeros(iN,1);

%% Matrix assembly
nElem = size(Elements,1);
iecount = 0;
for iElem = 1 : nElem
    elem = Elements(iElem,:);
    switch elem(2)
        %Normal finite elements
        case {12,23,24,34,36,38}
            dim = floor(elem(2)/10);                    % dimension number
            nodes = elem(4+(1:nNode(iElem)));           % elem nodes
            vert = model.Nodes(nodes,1+(1:dim)).';      % corner vertices
            [m, k] = elem_mk(elem(2), vert, model.Materials(elem(3),3:4));
            i = repmat(nodes.', 1, nNode(iElem));       % row indices
            j = i.';                                    % column indices
            ind = start(iElem)+1:start(iElem+1);
            I(ind) = i(:);
            J(ind) = j(:);
            M(ind) = m(:);
            K(ind) = k(:);
        %Infininite element
        case {111,122,133,134}
            %TODO: implementation should be nicer
            iecount = iecount +1;
            type = elem(2)-100;
            dim = floor(type/10);
            np = mod(type,10);%number of points at inner edge
            %Assemble matrix of mapping coordinates
            map = model.Nodes(elem(5:5+2*np-1),2:2+dim-1);
            P = model.Properties(elem(4),3);
            nm = (P-1)*np; %number of additional mapping nodes ( = (P-1)*np
            nodes = elem([5:5+np-1,5+2*np:5+2*np+nm-1]); %pressure node indices
            ptype = model.Properties(elem(4),4);
            param = model.Properties(elem(4),5:6);
            [m, k, c] = ielem_mkc(type,map,model.Materials(elem(3),3:4),P,ptype,param);
            i = repmat(nodes.', 1, inNode(iecount));       % row indices
            j = i.';      % column indices
            ind = istart(iecount)+1:istart(iecount+1);
            iI(ind) = i(:);
            iJ(ind) = j(:);
            iM(ind) = m(:);
            iK(ind) = k(:);
            iC(ind) = c(:);
    end
    progbar(1, nElem, iElem);
end
% Generate the sparse matrices from the arrays
nNodes = size(model.Nodes,1);
m = ismember(1:nNodes,[I;iI]);
vNodes =sum(m);
index = zeros(vNodes,1);
index(m) = 1:vNodes;
I = index(I);
J = index(J);
iI = index(iI);
iJ = index(iJ);
M = sparse(I, J, M, vNodes, vNodes);
K = sparse(I, J, K, vNodes, vNodes);

iM = sparse(iI, iJ, iM, vNodes, vNodes);
iK = sparse(iI, iJ, iK, vNodes, vNodes);
iC = sparse(iI, iJ, iC, vNodes, vNodes);

%% Postprocessing
% Ensure matrix symmetry
K = (K + K.')/2+iK;
M = (M + M.')/2+iM;
C = iC;
DOF = model.Nodes(m,1);
end