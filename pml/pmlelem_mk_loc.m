function [M_pml K_pml] = pmlelem_mk_loc(dim, xnode, w, N, dN, xi, om, hsigma, property, material)

ksi = (xi+1)/(2*property(2))+(property(3)-1)/property(2);   % ksi transformed to (0,1)
gamma = 1-1i/om*hsigma(ksi, material(2), property(1));      % gamma

ngauss = length(w);
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
M_pml = N.'*diag(q.*gamma)*N * material(1);
q = reshape(repmat((q./gamma).',dim,1),dim*ngauss,1);
K_pml = dNx.'*diag(q)*dNx * material(1)*material(2)^2;

end