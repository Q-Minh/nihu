function [M_pml K_pml] = pmlelem_mk_glob(type, xnode, w, N, dN, om, hsigma, psigma, L_pml, material)

dim = floor((type-200)/10);         % dimension number
nnode = mod((type-200),10);         % number of nodes of an element

% xyz(k,e,z)
xyz = N*xnode;
ngauss = length(w);

% gammax
if min(xyz(:,1)) < psigma{1}(1)
    gammax = 1-1i/om*hsigma{1}(abs(xyz(:,1)), abs(psigma{1}(1)-L_pml), material(2));
elseif max(xyz(:,1)) > psigma{1}(2)
    gammax = 1-1i/om*hsigma{1}(abs(xyz(:,1)), abs(psigma{1}(2)+L_pml), material(2));
else
    gammax = ones(ngauss,1);
end

% gammay
if isempty(hsigma{2})
    gammay = ones(ngauss,1);
elseif min(xyz(:,2)) < psigma{2}(1)
    gammay = 1-1i/om*hsigma{2}(abs(xyz(:,2)), abs(psigma{2}(1)-L_pml), material(2));
elseif max(xyz(:,2)) > psigma{2}(2)
    gammay = 1-1i/om*hsigma{2}(abs(xyz(:,2)), abs(psigma{2}(2)+L_pml), material(2));
else
    gammay = ones(ngauss,1);
end

% gammaz
if isempty(hsigma{3})
    gammaz = ones(ngauss,1);
elseif min(xyz(:,3)) < psigma{3}(1)
    gammaz = 1-1i/om*hsigma{3}(abs(xyz(:,3)), abs(psigma{3}(1)-L_pml), material(2));
elseif max(xyz(:,3)) > psigma{3}(2)
    gammaz = 1-1i/om*hsigma{3}(abs(xyz(:,3)), abs(psigma{3}(2)+L_pml), material(2));
else
    gammaz = ones(ngauss,1);
end

J = dN*xnode;          % Jacobian matrix
dNx = zeros(size(dN)); % dN(x,y,z)
j = zeros(ngauss,1);
for n = 1:ngauss
    ind = dim*(n-1)+(1:dim);
    Jac = J(ind,:);
    dNx(ind,:) = Jac\dN(ind,:);
    j(n) = abs(det(Jac));
end

% Numerical integration with Gaussian quadrature
q = w.*j;
M_pml = (repmat(gammax.*gammay.*gammaz,1,nnode).*N).'*diag(q)*N * material(1);
q = reshape(repmat(q.',dim,1),dim*ngauss,1);
g(:,:,1) = repmat(gammay.*gammaz./gammax,1,nnode);
if dim > 1
    g(:,:,2) = repmat(gammax.*gammaz./gammay,1,nnode);
end
if dim > 2
    g(:,:,3) = repmat(gammax.*gammay./gammaz,1,nnode);
end
K_pml = (reshape(shiftdim(g,2),numel(g)/nnode,nnode).*dNx).'*diag(q)*dNx * material(1)*material(2)^2;

end