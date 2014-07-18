function [N, dN, D, dD, dL, dMu] = ieshape34m(xi,P,ptype,a,b)
%JSHAPE2M calculates Jacobian shape functions of an element with radial
%order of Pand their gradients at given xi coordinates.

%TODO: Check for argument numbers
switch ptype
    case{'l','lag','lagrange'}
        z = linspace(-1,1-2/P,P);
        Lag = lagpolys(z);
        dLag = cell(P,1);
        % Scale the shape functions
        for c1 = 1:P;
            Lag{c1} = Lag{c1}/polyval(Lag{c1},z(c1))/((1-z(c1))/2);
            temp = Lag{c1}.*(length(Lag{c1})-1:-1:0);
            dLag{c1} = temp(1:end-1);
        end
    case{'j','jac','jacobi'}
    otherwise
        error('ieshape33m:invalidArg','Unknown polynomial family.');
end

s = xi(1,:).';
t = xi(2,:).';
u = xi(3,:).';
nxi = size(xi,2);

N=zeros(nxi,4*P);
dN = zeros(3*nxi,4*P);
dD = zeros(3*nxi,4*P);
dL = zeros(3*nxi,8);
dMu = zeros(3*nxi,4*P,4);

switch ptype
    %Lagrange shape functions
    case {'l','lag','lagrange'}
        for c1 = 1:P
            % Shape functions
            N(:,4*c1-3) = 1/8*(1-s).*polyval(Lag{c1},s).*(1-t).*(1-u);
            N(:,4*c1-2) = 1/8*(1-s).*polyval(Lag{c1},s).*(1+t).*(1-u);
            N(:,4*c1-1) = 1/8*(1-s).*polyval(Lag{c1},s).*(1+t).*(1+u);
            N(:,4*c1)   = 1/8*(1-s).*polyval(Lag{c1},s).*(1-t).*(1+u);
            % Derivatives
            % d/ds
            dN(1:3:end,4*c1-3) = 1/8*(1-t).*(1-u).*((-1)*polyval(Lag{c1},s)+(1-s).*polyval(dLag{c1},s));
            dN(1:3:end,4*c1-2) = 1/8*(1+t).*(1-u).*((-1)*polyval(Lag{c1},s)+(1-s).*polyval(dLag{c1},s));
            dN(1:3:end,4*c1-1) = 1/8*(1+t).*(1+u).*((-1)*polyval(Lag{c1},s)+(1-s).*polyval(dLag{c1},s));
            dN(1:3:end,4*c1)   = 1/8*(1-t).*(1+u).*((-1)*polyval(Lag{c1},s)+(1-s).*polyval(dLag{c1},s));
            % d/dt
            dN(2:3:end,4*c1-3) = 1/8*(1-s).*polyval(Lag{c1},s).*(1-u)*(-1);
            dN(2:3:end,4*c1-2) = 1/8*(1-s).*polyval(Lag{c1},s).*(1-u)*( 1);
            dN(2:3:end,4*c1-1) = 1/8*(1-s).*polyval(Lag{c1},s).*(1+u)*( 1);
            dN(2:3:end,4*c1)   = 1/8*(1-s).*polyval(Lag{c1},s).*(1+u)*(-1);
            % d/du
            dN(3:3:end,4*c1-3) = 1/8*(1-s).*polyval(Lag{c1},s).*(1-t)*(-1);
            dN(3:3:end,4*c1-2) = 1/8*(1-s).*polyval(Lag{c1},s).*(1+t)*(-1);
            dN(3:3:end,4*c1-1) = 1/8*(1-s).*polyval(Lag{c1},s).*(1+t)*( 1);
            dN(3:3:end,4*c1)   = 1/8*(1-s).*polyval(Lag{c1},s).*(1-t)*( 1);
        end

    %Jacobi shape functions
    case {'j','jac','jacobi'}
        for c1 = 1:P
            % Shape functions
            N(:,4*c1-3) = 1/8*(1-s).*jacpoly(c1-1,a,b,s).*(1-t).*(1-u);
            N(:,4*c1-2) = 1/8*(1-s).*jacpoly(c1-1,a,b,s).*(1+t).*(1-u);
            N(:,4*c1-1) = 1/8*(1-s).*jacpoly(c1-1,a,b,s).*(1+t).*(1+u);
            N(:,4*c1)   = 1/8*(1-s).*jacpoly(c1-1,a,b,s).*(1-t).*(1+u);
            % Derivatives
            % d/ds
            dN(1:3:end,4*c1-3) = 1/8*(1-t).*(1-u).*((-1)*jacpoly(c1-1,a,b,s)+(1-s).*djacpoly(c1-1,a,b,s,1));
            dN(1:3:end,4*c1-2) = 1/8*(1+t).*(1-u).*((-1)*jacpoly(c1-1,a,b,s)+(1-s).*djacpoly(c1-1,a,b,s,1));
            dN(1:3:end,4*c1-1) = 1/8*(1+t).*(1+u).*((-1)*jacpoly(c1-1,a,b,s)+(1-s).*djacpoly(c1-1,a,b,s,1));
            dN(1:3:end,4*c1)   = 1/8*(1-t).*(1+u).*((-1)*jacpoly(c1-1,a,b,s)+(1-s).*djacpoly(c1-1,a,b,s,1));
            % d/dt
            dN(2:3:end,4*c1-3) = 1/8*(1-s).*jacpoly(c1-1,a,b,s).*(1-u)*(-1);
            dN(2:3:end,4*c1-2) = 1/8*(1-s).*jacpoly(c1-1,a,b,s).*(1-u)*( 1);
            dN(2:3:end,4*c1-1) = 1/8*(1-s).*jacpoly(c1-1,a,b,s).*(1+u)*( 1);
            dN(2:3:end,4*c1)   = 1/8*(1-s).*jacpoly(c1-1,a,b,s).*(1+u)*(-1);
            % d/du
            dN(3:3:end,4*c1-3) = 1/8*(1-s).*jacpoly(c1-1,a,b,s).*(1-t)*(-1);
            dN(3:3:end,4*c1-2) = 1/8*(1-s).*jacpoly(c1-1,a,b,s).*(1+t)*(-1);
            dN(3:3:end,4*c1-1) = 1/8*(1-s).*jacpoly(c1-1,a,b,s).*(1+t)*( 1);
            dN(3:3:end,4*c1)   = 1/8*(1-s).*jacpoly(c1-1,a,b,s).*(1-t)*( 1);
        end
end
%Weighting function
D = repmat((1-s).^2/4,1,4*P);
%Gradient of the weighting function
dD(1:3:end,:) = repmat(-1/2+1/2*s,1,4*P);
dD(2:3:end,:) = zeros(nxi,4*P);
dD(3:3:end,:) = zeros(nxi,4*P);

%Gradients of mapping functions
dL(1:3:end,1) = -1/2*(1-u).*(1-t)./(s-1).^2;
dL(2:3:end,1) = -1/2*(1-u).*s./(s-1);
dL(3:3:end,1) = -1/2*(1-t).*s./(s-1);
dL(1:3:end,2) = -1/2*(1-u).*(1+t)./(s-1).^2;
dL(2:3:end,2) =  1/2*(1-u).*s./(s-1);
dL(3:3:end,2) = -1/2*(1+t).*s./(s-1);
dL(1:3:end,3) = -1/2*(1+u).*(1+t)./(s-1).^2;
dL(2:3:end,3) =  1/2*(1+u).*s./(s-1);
dL(3:3:end,3) =  1/2*(1+t).*s./(s-1);
dL(1:3:end,4) = -1/2*(1+u).*(1-t)./(s-1).^2;
dL(2:3:end,4) = -1/2*(1+u).*s./(s-1);
dL(3:3:end,4) =  1/2*(1-t).*s./(s-1);

dL(1:3:end,5) =  1/2*(1-u).*(1-t)./(s-1).^2;
dL(2:3:end,5) =  1/2*(1-u).*(1+s)./(2*(s-1));
dL(3:3:end,5) =  1/2*(1-t).*(1+s)./(2*(s-1));
dL(1:3:end,6) =  1/2*(1-u).*(1+t)./(s-1).^2;
dL(2:3:end,6) = -1/2*(1-u).*(1+s)./(2*(s-1));
dL(3:3:end,6) =  1/2*(1+t).*(1+s)./(2*(s-1));
dL(1:3:end,7) =  1/2*(1+u).*(1+t)./(s-1).^2;
dL(2:3:end,7) = -1/2*(1+u).*(1+s)./(2*(s-1));
dL(3:3:end,7) = -1/2*(1+t).*(1+s)./(2*(s-1));
dL(1:3:end,8) =  1/2*(1+u).*(1-t)./(s-1).^2;
dL(2:3:end,8) =  1/2*(1+u).*(1+s)./(2*(s-1));
dL(3:3:end,8) = -1/2*(1-t).*(1+s)./(2*(s-1));

%Gradients of phase functions
dMu(1:3:end,:,1) = repmat(1/2*(1-u-t+u.*t)./(s-1).^2,1,4*P);
dMu(2:3:end,:,1) = repmat(-1/4*(-1+u).*(1+s)./(s-1),1,4*P);
dMu(3:3:end,:,1) = repmat(-1/4*(-1+t).*(1+s)./(s-1),1,4*P);
dMu(1:3:end,:,2) = repmat(1/2*(1-u+t-u.*t)./(s-1).^2,1,4*P);
dMu(2:3:end,:,2) = repmat(-1/4*(1-u).*(1+s)./(s-1),1,4*P);
dMu(3:3:end,:,2) = repmat(-1/4*(-1-t).*(1+s)./(s-1),1,4*P);
dMu(1:3:end,:,3) = repmat(1/2*(1+u+t+u.*t)./(s-1).^2,1,4*P);
dMu(2:3:end,:,3) = repmat(-1/4*(1+u).*(1+s)./(s-1),1,4*P);
dMu(3:3:end,:,3) = repmat(-1/4*(1+t).*(1+s)./(s-1),1,4*P);
dMu(1:3:end,:,4) = repmat(1/2*(1+u-t-u.*t)./(s-1).^2,1,4*P);
dMu(2:3:end,:,4) = repmat(-1/4*(-1-u).*(1+s)./(s-1),1,4*P);
dMu(3:3:end,:,4) = repmat(-1/4*(1-t).*(1+s)./(s-1),1,4*P);