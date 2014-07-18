function [N, Mu, dN, dMu, dL, D, dD] = ieshape22m(xi,P,ptype,a,b)
%IESHAPE22ML calculates Lagrange shape functions of an element with radial
%order of P and their gradients at given xi coordinates.

error(nargchk(3,5,nargin,'struct'));
error(nargoutchk(1,7,nargout,'struct'));

%TODO: Check for argument numbers
switch ptype
    case{'l','lag','lagrange'}
        z = linspace(-1,1-2/P,P);
        Lag = lagpolys(z);
        dLag = cell(P,1);
        % Scale the shape functions
        for c1 = 1:P;
            Lag{c1} = Lag{c1}/polyval(Lag{c1},z(c1))/sqrt((1-z(c1))/2);
            temp = Lag{c1}.*(length(Lag{c1})-1:-1:0);
            dLag{c1} = temp(1:end-1);
        end
    case{'j','jac','jacobi'}
    otherwise
        error('ieshape22m:invalidArg','Unknown polynomial family.');
end

s = xi(:,1);
t = xi(:,2);
nxi = size(xi,1);

N=zeros(nxi,2*P);
if nargout > 2
    dN = zeros(2*nxi,2*P);
end

switch ptype
    %Lagrange polynomials
    case{'l','lag','lagrange'}
        for c1 = 1:P
            %Shape functions
            N(:,2*c1-1) = 1/(2*sqrt(2)).*sqrt(1-s).*polyval(Lag{c1},s).*(1-t);
            N(:,2*c1)   = 1/(2*sqrt(2)).*sqrt(1-s).*polyval(Lag{c1},s).*(1+t);
            %And their derivatives
            if nargout > 2
                %d/ds
                dN(1:2:end,2*c1-1) = -1/8*2^(1/2)./(1-s).^(1/2).*(1-t).*polyval(Lag{c1},s)+1/(2*sqrt(2)).*sqrt(1-s).*(1-t).*polyval(dLag{c1},s);
                dN(1:2:end,2*c1)   = -1/8*2^(1/2)./(1-s).^(1/2).*(1+t).*polyval(Lag{c1},s)+1/(2*sqrt(2)).*sqrt(1-s).*(1+t).*polyval(dLag{c1},s);
                %d/dt
                dN(2:2:end,2*c1-1) =  1/(2*sqrt(2)).*sqrt(1-s).*polyval(Lag{c1},s)*(-1);
                dN(2:2:end,2*c1)   =  1/(2*sqrt(2)).*sqrt(1-s).*polyval(Lag{c1},s)*( 1);
            end
        end
        %Jacobi polynomials
    case{'j','jac','jacobi'}
        for c1 = 1:P
            %Shape function
            N(:,2*c1-1) = 1/(2*sqrt(2)).*sqrt(1-s).*jacpoly(c1-1,a,b,s).*(1-t);
            N(:,2*c1)   = 1/(2*sqrt(2)).*sqrt(1-s).*jacpoly(c1-1,a,b,s).*(1+t);
            if nargout > 2
                %And their derivatives
                %d/ds
                dN(1:2:end,2*c1-1) = -1/8*2^(1/2)./(1-s).^(1/2).*(1-t).*jacpoly(c1-1,a,b,s)+1/(2*sqrt(2)).*sqrt(1-s).*(1-t).*djacpoly(c1-1,a,b,s,1);
                dN(1:2:end,2*c1)   = -1/8*2^(1/2)./(1-s).^(1/2).*(1+t).*jacpoly(c1-1,a,b,s)+1/(2*sqrt(2)).*sqrt(1-s).*(1+t).*djacpoly(c1-1,a,b,s,1);
                %d/dt
                dN(2:2:end,2*c1-1) =  1/(2*sqrt(2)).*sqrt(1-s).*jacpoly(c1-1,a,b,s)*(-1);
                dN(2:2:end,2*c1)   =  1/(2*sqrt(2)).*sqrt(1-s).*jacpoly(c1-1,a,b,s)*( 1);
            end
        end
end

if nargout > 1
    %Phase shape functions
    Mu = zeros(nxi,2);
    Mu(:,1) = 1/2*(1-t).*(1+s)./(1-s);
    Mu(:,2) = 1/2*(1+t).*(1+s)./(1-s);
end

if nargout > 3
    %Gradients of phase functions
    dMu = zeros(2*nxi,2);
    dMu(1:2:end,1) = (1/2-1/2*t)./(1-s)+(1/2-1/2*t).*(1+s)./(1-s).^2;
    dMu(2:2:end,1) = -1/2*(1+s)./(1-s);
    dMu(1:2:end,2) = (1/2+1/2*t)./(1-s)+(1/2+1/2*t).*(1+s)./(1-s).^2;
    dMu(2:2:end,2) = 1/2*(1+s)./(1-s);
end

if nargout > 4
    %Gradients of mapping functions
    dL = zeros(2*nxi,4);
    dL(1:2:end,1) = -(1-t)./(s-1).^2;
    dL(2:2:end,1) = -s./(s-1);
    dL(1:2:end,2) = -(1+t)./(s-1).^2;
    dL(2:2:end,2) = s./(s-1);
    dL(1:2:end,3) = (1-t)./(s-1).^2;
    dL(2:2:end,3) = (1+s)./(2*(s-1));
    dL(1:2:end,4) = (1+t)./(s-1).^2;
    dL(2:2:end,4) = -(1+s)./(2*(s-1));
end

if nargout > 5
    %Weighting function
    D = repmat((1-s).^2/4,1,2*P);
end

if nargout > 6
    %Gradient of the weighting function
    dD = zeros(2*nxi,2*P);
    dD(1:2:end,:) = repmat(-1/2+1/2*s,1,2*P);
    dD(2:2:end,:) = zeros(nxi,2*P);
end

end
