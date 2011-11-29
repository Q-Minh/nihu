P = 22;
alpha = 2;
beta = 0;
[ps, ws] = gaussquad(P+2);  % Gaussian quadrature in s direction
[ptu, wtu] = dunavant_rule(2);    % Gaussian quadrature in t direction
%Correct dunavant rules to match parent domain
ptu(1,:) = -2*ptu(1,:)+1;
ptu(2,:) =  2*ptu(2,:)-1;
ntu = size(ptu,2);
wtu      =  2*wtu;
ps = reshape(repmat(ps.',ntu,1),(P+2)*ntu,1);
ws = reshape(repmat(ws.',ntu,1),(P+2)*ntu,1);
ptu = repmat(ptu.',P+2,1);
wtu = repmat(wtu.',P+2,1);
xi = [ps(:) ptu(:,1) ptu(:,2)].';    % Gauss points as a list
w = ws.*wtu;                
%%
tic;
[N,dN,D,dD,dL,dMu] = ieshape33mj(xi,P,alpha,beta);
toc;
Map = [10 -1 -1; 10 1 -1; 10 1 1;...
       20 -1 -1; 20 1 -1; 20 1 1];
%%
tic;
[Ke,Ce,Me] = iemat33mj(Map,P,[P+2 2 2],alpha,beta);
toc;
tic;
a1 = norm(Map(1,:)-Map(4,:));
a2 = norm(Map(2,:)-Map(5,:));
a3 = norm(Map(3,:)-Map(6,:));
%% Secondary variables for coordinate transform
J = dL*Map;                     % Jacobian of the transformation
dMu = a1*dMu(:,:,1)+a2*dMu(:,:,2)+a3*dMu(:,:,3);          % Gradient of phse functions in stu space
dNx = zeros(size(dN));          % Gradient of shape functions in xyz space
dDx = zeros(size(dD));          % Gradient of weightinh functions in xyz space
dMux = zeros(size(dMu));        % Gradient of phase functions in xyz space
dMuxN = zeros(size(dMu));       % Gradient of phase functions x N
dDxN = zeros(size(dD));
dNxD = zeros(size(dD));
dMu2N = zeros(size(N));
dNdMu = zeros(size(N));
%% Coordinate transform
dim = 3;
ngauss = length(w);
j = zeros(ngauss,1);
for n = 1 : ngauss
    ind = dim*(n-1)+(1:dim);            %Current indices
    Jac = J(ind,:);                     %Current jacobian matrix
    dNx(ind,:)  = (Jac) \ dN(ind,:);
    dDx(ind,:)  = (Jac) \ dD(ind,:);
    dMux(ind,:) = (Jac) \ dMu(ind,:);
    diagN = diag(N(n,:));
    dMu2N(n,:)   = dot(dMux(ind,:),dMux(ind,:),1)*diagN;
    dNdMu(n,:)   = dot(dNx(ind,:),dMux(ind,:),1);
    dDxN(ind,:)  = dDx(ind,:)*diagN;
    dMuxN(ind,:) = dMux(ind,:)*diagN;
    dNxD(ind,:)  = dNx(ind,:)*diag(D(n,:));
    j(n) = abs(det(Jac));
end
%% Numerical integration with Gaussian quadrature
q = w.*j;
qd = reshape(repmat(q.',dim, 1),dim*ngauss,1);
K = (dNx.'* diag(qd)*dDxN+dNx.'*diag(qd) *dNxD).';
C = (dNdMu.'*diag(q)*(D.*N)).'-(dMuxN.'*diag(qd)*dDxN).'-(dMuxN.'*diag(qd)*dNxD).';
M = N.'*diag(q)*(D.*N) - dMu2N.'*diag(q)*(D.*N);
toc;
%%
disp('Errors:')
disp(['K: ',num2str(norm(K-Ke))]);
disp(['C: ',num2str(norm(C-Ce))]);
disp(['M: ',num2str(norm(M-Me))]);