function [A DOF] = model_a_cyl(model)
%MODEL_A_CYL Compute excitation surface matrix of the FE model
%   [A DOF] = MODEL_A_CYL(MODEL) Computes the boundary excitation matrix
%   A of the cylindrically symmetric finite element model MODEL.
%
% See also: elem_mk, elem_mk_cyl, elem_a, elem_a_cyl, model_mk,
% model_mk_cyl, model_a

% Peter Fiala
% 2008 May

%% Preconditioning of the model
Elements = drop_IDs(model);
nElem = size(Elements,1);
Nodes = model.Nodes;
nNode = size(Nodes,1);

%% Matrix assembly
A = sparse(nNode, nNode);
wb = waitbar(0, 'Assembling model excitaion matrix...');
for iElem = 1 : nElem
    row = Elements(iElem,:);
    nod = row(5:6);
    xnode = Nodes(nod,2:4).';
    a = elem_a_cyl(xnode, model.Materials(row(3),[3 4]));
    A(nod,nod) = A(nod,nod) + a;
    waitbar(iElem/nElem, wb);
end
close(wb);

%% Postprocessing
% Ensure matrix symmetry
A = (A + A.')/2;
% Degree of Freedom vector
DOF = Nodes(:,1);

end

function [A] = elem_a_cyl(xnode, material)
%ELEM_A_CYL  Acoustic excitation matrices for cylindrical geometries
%   A = ELEM_A_CYL(XNODE, MAT) Computes the element excitation
%   matrix A of an acoustic 2-node LINE elem.
%   XNODE : 2x2 matrix containing the rz coordinates of each node
%   MAT   : [rho c] vector containing the mass density and sound volocity
%   A     : 2x2 element excitation matrix
%
% See also: elem_mk, elem_mk_cyl, elem_a, model_mk, model_mk_cyl,
% model_a, model_a_cyl

% Fiala Peter
% 2008. may

dr = diff(xnode,1,2);
r = xnode(1,:);
A = (diag(r) + mean(r))*sqrt(dot(dr,dr,1))/6*(material(1)^2*material(2)^2);
end