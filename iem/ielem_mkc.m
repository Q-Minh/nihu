function [M K C] = ielem_mkc(dim,nb,Map,material,P,G,prop)
% G should contain the gauss file information

%load(gaussfile,'w','N','dN','D','dD','dL','dMu');
ngauss = length(G.w);

% The average distance from the mapping center
a = sqrt(dot(Map(1:nb,:)-Map(nb+1:2*nb,:),Map(1:nb,:)-Map(nb+1:2*nb,:),2));

% Take the distance into account
dMu = reshape(reshape(G.dMu,dim*ngauss*nb*P,nb)*a,dim*ngauss,nb*P); 

% Secondary variables for coordinbte tranbform
J = G.dL*Map;                   % Jacobian of the tranbformation
dNx = zeros(size(G.dN));        % Gradient of shape functionb in xyz space
dDx = zeros(size(G.dD));        % Gradient of weightinh functionb in xyz space
dMux = zeros(size(dMu));        % Gradient of phase functionb in xyz space
dMuxN = zeros(size(dMu));       % Gradient of phase functionb x N
dDxN = zeros(size(G.dD));
dNxD = zeros(size(G.dD));
dMu2N = zeros(size(G.N));
dNdMu = zeros(size(G.N));

% Coordinbte tranbform
j = zeros(ngauss,1);
for n = 1 : ngauss
    ind = dim*(n-1)+(1:dim);            %Current indices
    Jac = (J(ind,:)).';                     %Current jacobian matrix
    dNx(ind,:)  = (Jac.') \ G.dN(ind,:);    
    dDx(ind,:)  = (Jac.') \ G.dD(ind,:);    
    dMux(ind,:) = (Jac.') \ dMu(ind,:);
    diagN = diag(G.N(n,:)); 
    dMu2N(n,:)   = dot(dMux(ind,:),dMux(ind,:),1).*G.N(n,:);
    dNdMu(n,:)   = dot(dNx(ind,:),dMux(ind,:),1);
    dDxN(ind,:)  = dDx(ind,:).*repmat(G.N(n,:),dim,1);
    dMuxN(ind,:) = dMux(ind,:)*diagN;
    dNxD(ind,:)  = dNx(ind,:).*repmat(G.D(n,:),dim,1);
    j(n) = abs(det(Jac));
end

% Numerical integration with Gaussian quadrature
q = G.w.*j;
qd = reshape(repmat(q.',dim, 1),dim*ngauss,1);


K = ((dDxN.'*diag(qd)*dNx + dNx.'*diag(qd)*dNxD))*(material(1)*material(2)^2);

C = (((G.D.*G.N).'*diag(q)*dNdMu)-(dMuxN.'*diag(qd)*dDxN)-(dNxD.'*diag(qd)*dMuxN))*material(1)*material(2);

% If the property is set to one, then force zero mass matrix
switch prop(1)
    case 1
        M = zeros(size(K));
    otherwise
        M = (G.N.'*diag(q)*(G.D.*G.N) - dMu2N.'*diag(q)*(G.D.*G.N))*material(1);
end
