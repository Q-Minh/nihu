function [N, Mu, dN, dMu, dL, D, dD] = ie_shapefun(type,xi,P,varargin)
%IESHAPE22ML calculates Lagrange shape functions of an element with radial
%order of P and their gradients at given xi coordinates.

% TODO: 1D linear infinite elements (type == 111)

% Check arguments
error(nargchk(3,6,nargin,'struct'));
error(nargoutchk(1,7,nargout,'struct'));

% Check element type
switch type
    case 111
        dim = 1;
        s = xi(:,1);
    case 122
        dim = 2;
        s = xi(:,1);
        t = xi(:,2);
    case {133, 134}
        dim = 3;
        s = xi(:,1);
        t = xi(:,2);
        u = xi(:,3);
    otherwise
        error('ie_shapefun:invalidArg',...
            ['Unknown inifinite element type: ',num2str(type),'.']);
end

% Number of base points
base = mod(type,10);

% If polynomial type not defined use Lagrange
if numel(varargin) == 0
    ptype = 'l';
else
    ptype = varargin{1};
end

% Pre-create polynomials for shape functions
switch ptype
    case{'l','lag','lagrange'}
        z = linspace(-1,1-2/P,P);
        Lag = lagpolys(z);
        dLag = cell(P,1);
        % Scale the shape functions
        for c1 = 1:P;
            switch type
                case 122
                    Lag{c1} = Lag{c1}/polyval(Lag{c1},z(c1))/sqrt((1-z(c1))/2);
                case {133, 134}
                    Lag{c1} = Lag{c1}/polyval(Lag{c1},z(c1))/((1-z(c1))/2);
            end
            temp = Lag{c1}.*(length(Lag{c1})-1:-1:0);
            dLag{c1} = temp(1:end-1);
        end
    case{'j','jac','jacobi'}
        a = varargin{2};
        b = varargin{3};
    otherwise
        error('ie_shapefun:invalidArg','Unknown polynomial family.');
end

% Number of xi points
nxi = size(xi,1);

% Pressure shape functions
N=zeros(nxi,base*P);
if nargout > 2
    dN = zeros(dim*nxi,base*P);
end

switch type
    % 2D quad infinite elements
    case 122
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
    case 133
        % 3D penta infinite elements
        switch ptype
            %Lagrange shape functions
            case {'l','lag','lagrange'}
                for c1 = 1:P
                    % Shape functions
                    N(:,3*c1-2) = 1/4*(1-s).*polyval(Lag{c1},s).*(1-t);
                    N(:,3*c1-1) = 1/2*(1-s).*polyval(Lag{c1},s).*(1-(1+u)/2-(1-t)/2);
                    N(:,3*c1)   = 1/4*(1-s).*polyval(Lag{c1},s).*(1+u);
                    % Calculate gradients
                    if nargout > 2
                        % d/ds
                        dN(1:3:end,3*c1-2) =               1/4*(1-t).*((-1)*polyval(Lag{c1},s)+(1-s).*polyval(dLag{c1},s));
                        dN(1:3:end,3*c1-1) = 1/2*(1-(1+u)/2-(1-t)/2).*((-1)*polyval(Lag{c1},s)+(1-s).*polyval(dLag{c1},s));
                        dN(1:3:end,3*c1)   =               1/4*(1+u).*((-1)*polyval(Lag{c1},s)+(1-s).*polyval(dLag{c1},s));
                        % d/dt
                        dN(2:3:end,3*c1-2) = 1/4*(1-s).*polyval(Lag{c1},s)*(-1);
                        dN(2:3:end,3*c1-1) = 1/4*(1-s).*polyval(Lag{c1},s)*( 1);
                        dN(2:3:end,3*c1)   =                                   0;
                        % d/du
                        dN(3:3:end,3*c1-2) = 0;
                        dN(3:3:end,3*c1-1) = 1/4*(1-s).*polyval(Lag{c1},s)*(-1);
                        dN(3:3:end,3*c1)   = 1/4*(1-s).*polyval(Lag{c1},s)*( 1);
                    end
                end
                %Jacobi shape functions
            case {'j','jac','jacobi'}
                for c1 = 1:P
                    % Shape functions
                    N(:,3*c1-2) = 1/4*(1-s).*jacpoly(c1-1,a,b,s).*(1-t);
                    N(:,3*c1-1) = 1/2*(1-s).*jacpoly(c1-1,a,b,s).*(1-(1+u)/2-(1-t)/2);
                    N(:,3*c1)   = 1/4*(1-s).*jacpoly(c1-1,a,b,s).*(1+u);
                    % Calculate gradients
                    if nargout > 2
                        % d/ds
                        dN(1:3:end,3*c1-2) =               1/4*(1-t).*((-1)*jacpoly(c1-1,a,b,s)+(1-s).*djacpoly(c1-1,a,b,s,1));
                        dN(1:3:end,3*c1-1) = 1/2*(1-(1+u)/2-(1-t)/2).*((-1)*jacpoly(c1-1,a,b,s)+(1-s).*djacpoly(c1-1,a,b,s,1));
                        dN(1:3:end,3*c1)   =               1/4*(1+u).*((-1)*jacpoly(c1-1,a,b,s)+(1-s).*djacpoly(c1-1,a,b,s,1));
                        % d/dt
                        dN(2:3:end,3*c1-2) = 1/4*(1-s).*jacpoly(c1-1,a,b,s)*(-1);
                        dN(2:3:end,3*c1-1) = 1/4*(1-s).*jacpoly(c1-1,a,b,s)*( 1);
                        dN(2:3:end,3*c1)   =                                   0;
                        % d/du
                        dN(3:3:end,3*c1-2) = 0;
                        dN(3:3:end,3*c1-1) = 1/4*(1-s).*jacpoly(c1-1,a,b,s)*(-1);
                        dN(3:3:end,3*c1)   = 1/4*(1-s).*jacpoly(c1-1,a,b,s)*( 1);
                    end
                end
        end
    case 134
        % 3D hexa infinite elements
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
                    if nargout > 2
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
                    if nargout > 2
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
        end
end

%Phase shape functions
Mu = zeros(nxi,base);
if nargout > 1
    switch type
        case 122
            Mu(:,1) = 1/2*(1-t).*(1+s)./(1-s);
            Mu(:,2) = 1/2*(1+t).*(1+s)./(1-s);
        case 133
            Mu(:,1) = 1/2*(1-t).*(1+s)./(1-s);
            Mu(:,2) = (1-(1+u)/2-(1-t)/2).*(1+s)./(1-s);
            Mu(:,3) = 1/2*(1+u).*(1+s)./(1-s);
        case 134
            Mu(:,1) = 1/4*(1-t).*(1-u).*(1+s)./(1-s);
            Mu(:,2) = 1/4*(1+t).*(1-u).*(1+s)./(1-s);
            Mu(:,3) = 1/4*(1+t).*(1+u).*(1+s)./(1-s);
            Mu(:,4) = 1/4*(1-t).*(1+u).*(1+s)./(1-s);
    end
end

%Gradients of phase functions
dMu = zeros(dim*nxi,base);
if nargout > 3
    switch type
        case 122
            dMu(1:2:end,1) = (1/2-1/2*t)./(1-s)+(1/2-1/2*t).*(1+s)./(1-s).^2;
            dMu(2:2:end,1) = -1/2*(1+s)./(1-s);
            dMu(1:2:end,2) = (1/2+1/2*t)./(1-s)+(1/2+1/2*t).*(1+s)./(1-s).^2;
            dMu(2:2:end,2) = 1/2*(1+s)./(1-s);
        case 133
            dMu(1:3:end,1) =     (1-t)./(s-1).^2;
            dMu(2:3:end,1) = 1/2*(1+s)./(s-1);
            dMu(3:3:end,1) = 0;
            dMu(1:3:end,2) = (1+u)./(s-1).^2;
            dMu(2:3:end,2) = 0;
            dMu(3:3:end,2) = -1/2*(1+s)./(s-1);
            dMu(1:3:end,3) = (-u+t)./(s-1).^2;
            dMu(2:3:end,3) = -1/2*(1+s)./(s-1);
            dMu(3:3:end,3) =  1/2*(1+s)./(s-1);
        case 134
            dMu(1:3:end,1) = 1/2*(1-u-t+u.*t)./(s-1).^2;
            dMu(2:3:end,1) = -1/4*(-1+u).*(1+s)./(s-1);
            dMu(3:3:end,1) = -1/4*(-1+t).*(1+s)./(s-1);
            dMu(1:3:end,2) = 1/2*(1-u+t-u.*t)./(s-1).^2;
            dMu(2:3:end,2) = -1/4*(1-u).*(1+s)./(s-1);
            dMu(3:3:end,2) = -1/4*(-1-t).*(1+s)./(s-1);
            dMu(1:3:end,3) = 1/2*(1+u+t+u.*t)./(s-1).^2;
            dMu(2:3:end,3) = -1/4*(1+u).*(1+s)./(s-1);
            dMu(3:3:end,3) = -1/4*(1+t).*(1+s)./(s-1);
            dMu(1:3:end,4) = 1/2*(1+u-t-u.*t)./(s-1).^2;
            dMu(2:3:end,4) = -1/4*(-1-u).*(1+s)./(s-1);
            dMu(3:3:end,4) = -1/4*(1-t).*(1+s)./(s-1);
    end
end

%Gradients of mapping functions
dL = zeros(dim*nxi,2*base);
if nargout > 4
    switch type
        % 2D quad infinite element
        case 122
            dL(1:2:end,1) = -(1-t)./(s-1).^2;
            dL(2:2:end,1) = -s./(s-1);
            dL(1:2:end,2) = -(1+t)./(s-1).^2;
            dL(2:2:end,2) = s./(s-1);
            dL(1:2:end,3) = (1-t)./(s-1).^2;
            dL(2:2:end,3) = (1+s)./(2*(s-1));
            dL(1:2:end,4) = (1+t)./(s-1).^2;
            dL(2:2:end,4) = -(1+s)./(2*(s-1));
        % 3D penta infinite element
        case 133
            dL(1:3:end,1) = -(1-t)./(s-1).^2;
            dL(2:3:end,1) = -s./(s-1);
            dL(3:3:end,1) =  0;
            dL(1:3:end,2) =  (u-t)./(s-1).^2;
            dL(2:3:end,2) =  s./(s-1);
            dL(3:3:end,2) = -s./(s-1);
            dL(1:3:end,3) = -(1+u)./(s-1).^2;
            dL(2:3:end,3) =  0;
            dL(3:3:end,3) =  s./(s-1);
            
            dL(1:3:end,4) = (1-t)./(s-1).^2;
            dL(2:3:end,4) = 1/2*(1+s)./(s-1);
            dL(3:3:end,4) = 0;
            dL(1:3:end,5) = -(u-t)./(s-1).^2;
            dL(2:3:end,5) = -1/2*(1+s)./(s-1);
            dL(3:3:end,5) =  1/2*(1+s)./(s-1);
            dL(1:3:end,6) = (1+u)./(s-1).^2;
            dL(2:3:end,6) = 0;
            dL(3:3:end,6) = -1/2*(1+s)./(s-1);
        % 3D hexa infinite element
        case 134
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
    end
end

%Weighting function
if nargout > 5
    D = repmat((1-s).^2/4,1,base*P);
end

%Gradient of the weighting function
if nargout > 6
    dD = zeros(dim*nxi,base*P);
    
    switch type
        case 122
            dD(1:2:end,:) = repmat(-1/2+1/2*s,1,base*P);
            dD(2:2:end,:) = zeros(nxi,base*P);
        case {133,134}
            dD(1:3:end,:) = repmat(-1/2+1/2*s,1,base*P);
            dD(2:3:end,:) = zeros(nxi,base*P);
            dD(3:3:end,:) = zeros(nxi,base*P);
    end
end

end
