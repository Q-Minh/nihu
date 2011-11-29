function [A DOF] = model_a(boundary)
%boundary_A Compute excitation surface matrix of the FE boundary
%   [A DOF] = model_A(boundary) Computes the boundary excitation matrix A of
%   the finite element boundary boundary.
%
% See also: elem_mk, elem_mk_cyl, elem_a, elem_a_cyl, model_mk,
% model_mk_cyl, model_a_cyl

% Copyright 2010 Peter Fiala
% 2009

%% Preprocessing
Elements = drop_IDs(boundary); % drop IDs

%% Preallocating space for the acoustic matrices
% The elements of the sparse matrices are collected in arrays. The matrix
% values of each element are stored in non overlapping parts of the arrays.
nNode = sum(Elements(:,5:end)~=0, 2);   % number of nodes per element
nVal = nNode.^2;                        % number of matrix entries per elem
start = [0; cumsum(nVal)];              % starting indices of elements
N = start(end); % total number of matrix entries (non overlapped)
I = zeros(N,1); % preallocating row indices
J = zeros(N,1); % preallocating column indices
A = zeros(N,1); % preallocating mass values

%% Matrix assembly
nElem = size(Elements,1);
for iElem = 1 : nElem
    elem = Elements(iElem,:);
    dim = floor(elem(2)/10)+1;                  % dimension number
    nodes = elem(4+(1:nNode(iElem)));           % elem nodes
    vert = boundary.Nodes(nodes,1+(1:dim)).';      % corner vertices
    a = elem_a(elem(2), vert, boundary.Materials(elem(3),[3 4]));
    i = repmat(nodes.', 1, nNode(iElem));       % row indices
    j = i.';                                    % column indices
    ind = start(iElem)+1:start(iElem+1);
    I(ind) = i(:);
    J(ind) = j(:);
    A(ind) = a(:);
    progbar(1, nElem, iElem);
end
nNodes = size(boundary.Nodes,1);
A = sparse(I, J, A, nNodes, nNodes);

%% Postprocessing
% Ensure matrix symmetry
A = (A + A.')/2;
% Degree of Freedom vector
DOF = boundary.Nodes(:,1);
end

function [A] = elem_a(type, xnode, material)
%ELEM_A  Acoustic excitation matrices
%   A = ELEM_A(TYPE, XNODE, MAT) Computes the element excitation
%   matrix A of an acoustic element.
%   TYPE  : Element type selector [12 23 24]
%   XNODE : nDim x nNodes matrix containing the coordinates of each node.
%   MAT   : [rho c] vector containing the mass density and sound velocity
%   A     : nNodes x nNodes excitation matrix
%
% See also: elem_mk, elem_mk_cyl, elem_a_cyl, model_mk, model_mk_cyl,
% model_a, model_a_cyl

% Fiala Peter
% 2008. may

switch type
    case 12
        %% line element (boundary of 2D model)
        dr = diff(xnode,1,2);
        j = sqrt(dot(dr,dr,1));
        A = [
            2 1
            1 2
            ] * (j/6*material(1)^2*material(2)^2);
    case 23
        %% tria element (boundary of 3D model)
        df = cross(diff(xnode(:,[1 2]), 1, 2), diff(xnode(:,[1 3]), 1, 2));
        A = [
            2 -1 -1
            -1 2 3
            -1 3 6
            ]*sqrt(dot(df,df))/24*(material(1)^2*material(2)^2);
    case 24
        %% quad element (boundary of 3D model)
        load(fullfile(fileparts(mfilename('fullpath')),...
            'data', 'gauss_quad_3'), 'w', 'N', 'dNxi');
        ngauss = length(w);
        J = xnode * dNxi.';
        df = zeros(3,ngauss);
        for n = 1 : ngauss
            df(:,n) = cross(J(:,2*n-1), J(:,2*n));
        end
        j = sqrt(dot(df,df,1)).';
        A = (N.'*diag(w.*j)*N)*(material(1)^2*material(2)^2);
end
end