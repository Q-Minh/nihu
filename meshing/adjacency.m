function [S, DOF] = adjacency(mesh)
%ADJACENCY  Compute adjacency matrix of a fe mesh
%   [S, DOF] = adjacency(mesh) computes the adjacency matrix of a mesh. The
%   adjacency matrix is based on the physical ordering of the nodes and not
%   on the node IDs.
% Parameters:
%   mesh
% Output:
%   S   : NxN sparse adjacency matrix
%   DOF : The node IDs corresponding to the matrix elements

%   Copyright 2008-2010 P. Fiala
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last modified: 2012.12.10.

%% Preprocessing
% drop User IDs
Elements = drop_IDs(mesh);

%% Preallocating space for the adjacency matrices
% The elements of the sparse matrices are collected in arrays. The matrix
% values of each element are stored in non overlapping parts of the arrays.
nNode = sum(Elements(:,5:end)~=0, 2);   % number of nodes per element
nVal = nNode.^2;                        % number of matrix entries per elem
start = [0; cumsum(nVal)];              % starting indices of elements
N = start(end); % total number of matrix entries (non overlapped)
I = zeros(N,1); % preallocating row indices
J = zeros(N,1); % preallocating column indices

%% Filling elements of S
nElem = size(Elements,1);
for iElem = 1 : nElem
    elem = Elements(iElem,:);
    nodes = elem(4+(1:nNode(iElem)));           % elem nodes
    i = repmat(nodes.', 1, nNode(iElem));       % row indices
    j = i.';                                    % column indices
    ind = start(iElem)+1:start(iElem+1);
    I(ind) = i(:);
    J(ind) = j(:);
    progbar(1, nElem, iElem);
end
% Generate the sparse matrices from the arrays
nNodes = size(mesh.Nodes,1);
S = sparse(I, J, ones(size(I)), nNodes, nNodes);
S = logical(S);

%% Postprocessing
% Degree of Freedom vector
DOF = mesh.Nodes(:,1);
end
