function [S, E] = adjacency(mesh)
%ADJACENCY  Compute adjacency matrix of a fe mesh
%   [S, E] = adjacency(mesh) computes the adjacency matrices of a mesh.
%   The adjacency matrices are based on the physical ordering of the nodes
%   and not on the node IDs.
% Parameters:
%   mesh
% Output:
%   S : nNode x nNode sparse adjacency matrix
%   E : nNode x nElem sparse connectivity matrix

%   Copyright 2008-2012 P. Fiala, P. Rucz
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last modified: 2012.12.19.

nNodes = size(mesh.Nodes,1);
nElem = size(mesh.Elements,1);

% drop User IDs
Elements = drop_IDs(mesh);

% Preallocating space for the adjacency matrices
% The elements of the sparse matrices are collected in arrays. The matrix
% values of each element are stored in non overlapping parts of the arrays.
nNode = sum(Elements(:,5:end)~=0, 2);   % number of nodes per element
nVal = nNode.^2;                        % number of matrix entries per elem
start = [0; cumsum(nVal)];              % starting indices of elements
N = start(end); % total number of matrix entries (non overlapped)
I = zeros(N,1); % preallocating row indices
J = zeros(N,1); % preallocating column indices

% Filling elements of S
for iElem = 1 : nElem
    elem = Elements(iElem,:);
    nodes = elem(4+(1:nNode(iElem)));           % elem nodes
    i = repmat(nodes.', 1, nNode(iElem));       % row indices
    j = i.';                                    % column indices
    ind = start(iElem)+1:start(iElem+1);
    I(ind) = i(:);
    J(ind) = j(:);
%     progbar(1, nElem, iElem);
end
% Generate the sparse matrices from the arrays
S = sparse(I, J, ones(size(I)), nNodes, nNodes);
S = logical(S);

if nargout > 1
    % Filling elements of E
    e = Elements(:,5:end);
    i = repmat((1:nElem)',1,size(e,2));
    v = ones(size(e));
    E = sparse(i(e~=0), e(e~=0), v(e~=0), nElem, nNodes);
    E = logical(E);
end

end
