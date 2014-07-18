function [Ap, Kp, Av, Kv] = ie_postproc(type,Map,P,xi)

% Check arguments
error(nargchk(4,4,nargin,'struct'));
error(nargoutchk(2,4,nargout,'struct'));

% Get element base
base = mod(type,10);

%Calculate distance of nodes from location of virtual source points
a = sqrt(dot(Map(1:base,:),Map(base+1:2*base,:),2));

% Create pressure interpolation function
% p_interp = (Ap.*exp(-1i*k*Kp))*p (be careful with exp(Kp)!)
if nargout < 3
  [Ap, Mu] = ie_shapefun(type,xi,P,'l');
  Kp = repmat(Mu*a,1,base*P);
else
    nxi = size(xi,1);
    dim = ceil(log2(base))+1;
    [Ap, Mu, dN, dMu, dL] = ie_shapefun(type,xi,P,'l');
    Kp = repmat(Mu*a,1,base*P);
    dMu = repmat(dMu*a,1,base*P);
    J = dL*Map;
    for c1 = 1:nxi
        Av(dim*(c1-1)+1:dim*c1,:) =  J(dim*(c1-1)+1:dim*c1,:) \ dN(dim*(c1-1)+1:dim*c1,:);
        Kv(dim*(c1-1)+1:dim*c1,:) =  J(dim*(c1-1)+1:dim*c1,:) \ dMu(dim*(c1-1)+1:dim*c1,:);
    end
end

end