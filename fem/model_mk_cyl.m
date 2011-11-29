function [M, K, B, DOF] = model_mk_cyl(model)
%MODEL_MK_CYL Compute mass and stiffness matrices of the FE model
%   [M K B DOF] = MODEL_MK_CYL(MODEL) Computes the mass matrix M and
%   stiffness matrices K and B of a FE model MODEL. The output vector DOF
%   relates the elements of the output matrices to the degrees of freedom
%   of the model.
%
% See also: elem_mk, elem_mk_cyl, elem_a, elem_a_cyl, model_mk,
% model_mk, model_a, model_a_cyl

% Peter Fiala
% 2008 May

%% Preconditioning
% Drop non-acoustical elements
AcMat = model.Materials(model.Materials(:,2) == 1,1);
model.Elements = model.Elements(ismember(model.Elements(:,3), AcMat),:);
% drop User IDs
model.Elements = drop_IDs(model);

%% Computing the acoustic matrices
nNode = size(model.Nodes,1);
nElem = size(model.Elements,1);
M = sparse(nNode,nNode);
K = sparse(nNode,nNode);
B = sparse(nNode,nNode);
wb = waitbar(0, 'Assembling model mass and stiffness matrices...');
for iElem = 1 : nElem
    elem = model.Elements(iElem,:);
    nNodes = mod(elem(2), 10);
    nodes = elem(4+(1:nNodes));
    xnode = model.Nodes(nodes,[2 4]).';
    [m, k, b] = elem_mk_cyl(elem(2), xnode, model.Materials(elem(3),3:4));
    M(nodes,nodes) = M(nodes,nodes) + m;
    K(nodes,nodes) = K(nodes,nodes) + k;
    B(nodes,nodes) = B(nodes,nodes) + b;
    waitbar(iElem/nElem, wb);
end
close(wb);

%% Postprocessing
% Ensure matrix symmetry
K = (K + K.')/2;
M = (M + M.')/2;
B = (B + B.')/2;
% r=0 DOF of B must vanish
ind = find(model.Nodes(:,2) < eps);
B(ind,ind) = 0;
% Degree of Freedom vector
DOF = model.Nodes(:,1);

end

function [M, K, B] = elem_mk_cyl(type, xnode, material)
%ELEM_MK_CYL  Acoustic element matrices for cylindrical geometries
%   [M, K] = ELEM_MK_ACOU_CYL_QUAD(TYPE, XNODE, MAT)
%   [M, K, B] = ELEM_MK_ACOU_CYL_QUAD(TYPE, XNODE, MAT)
%   Computes the element matrices of an acoustic elem.
%   TYPE  : Element type selector (23, 24)
%   XNODE : 2 x nNodes matrix containing the rz coordinates of each node
%   MAT   : [rho c] vector containing the mass density and sound velocity
%   M,K,B : nNodes x nNodes element mass and stiffness matrices
%
% See also: elem_mk, elem_a, elem_a_cyl, model_mk, model_mk_cyl,
% model_a, model_a_cyl

% Fiala Peter
% 2008. June

%% Load Gaussian quadrature weights and shape function samples
switch type
    case 23
        load data/gauss_tria_5 w N dNxi
    case 24
        load data/gauss_quad_5 w N dNxi
end
ngauss = length(w);

%% Coordinate transform
J = xnode * dNxi.';
dNx = zeros(size(dNxi));
r = N * xnode(1,:).';

j = zeros(ngauss,1);
for n = 1 : ngauss
    ind = 2*(n-1)+(1:2);
    Jac = J(:,ind);
    dNx(ind,:) = (Jac.') \ dNxi(ind,:);
    j(n) = abs(det(Jac));
end

%% Numerical integration
q1 = w.*j.*r;
M = N.'*diag(q1)*N * material(1);
q2 = w.*j./r;
B = N.'*diag(q2)*N * (material(1)*material(2)^2);
q1 = reshape(repmat(q1.',2, 1),2*ngauss,1);
K = dNx.' * diag(q1) * dNx * (material(1)*material(2)^2);

end