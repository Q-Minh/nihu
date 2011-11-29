function [Ke,Ce,Me] = iemat2m(M,z,ngauss)
%IEMAT2M Computes element matrices for 2D infinite element.
%   M      : Coordinates of the mapping nodes, must be a 4x2 matrix.
%   Z      : Zero points of the radial shape functions in the parent
%             coordinate system. The length of Z determines the radial
%             order of the element and therefore the size of element
%             matrices.
%   ngauss : Number of base points of gaussian quadrature integration.
%
%   Ke     : Element stiffness matrix.
%   Ce     : Element damping matrix.
%   Me     : Element mass matrix.
%
% See also: iemap2D, lagpolys, gaussquad, iepproc2m

% Peter Rucz
% 2009 november


%Create lagrangian shape functions 
l = length(z);
L = lagpolys(z);
dL = cell(l,1);
% Scale the shape functions
for c1 = 1:length(z);
    L{c1} = L{c1}/polyval(L{c1},z(c1))/sqrt(1-z(c1));
    temp = L{c1}.*(length(L{c1})-1:-1:0);
    dL{c1} = temp(1:end-1);
end

% Coordinates of the mapping
x1=M(1,1); y1=M(1,2);
x2=M(2,1); y2=M(2,2);
x3=M(3,1); y3=M(3,2);
x4=M(4,1); y4=M(4,2);

%Calculate distance of nodes from location of virtual source points
a1 = norm(M(1,:)-M(3,:));
a2 = norm(M(2,:)-M(4,:));

%Gaussian quadrature base points and weights
[ps,ws] = gaussquad(ngauss(1));
[pt,wt] = gaussquad(ngauss(2));

Ke = zeros(2*l);
Ce = zeros(2*l);
Me = zeros(2*l);

N = zeros(1,2*l);
gN = zeros(2,2*l);

for c1 = 1:ngauss(1)
    for c2 = 1:ngauss(2)
        s = ps(c1);
        t = pt(c2);
        %Calculate shape functions
        for c3 = 1:l
            N(2*(c3-1)+1) = 1/2*sqrt(1-s)*polyval(L{c3},s)*(1-t);
            N(2*(c3-1)+2) = 1/2*sqrt(1-s)*polyval(L{c3},s)*(1+t);
            gN(1,2*(c3-1)+1) = 1/2*sqrt(1-s)*(1-t)*polyval(dL{c3},s)-1/4/(1-s)^(1/2)*(1-t)*polyval(L{c3},s);
            gN(1,2*(c3-1)+2) = 1/2*sqrt(1-s)*(1+t)*polyval(dL{c3},s)-1/4/(1-s)^(1/2)*(1+t)*polyval(L{c3},s);
            gN(2,2*(c3-1)+1) = -1/2*sqrt(1-s)*polyval(L{c3},s);
            gN(2,2*(c3-1)+2) =  1/2*sqrt(1-s)*polyval(L{c3},s);
        end
        %Weighting function
        D = (1-s)^2/4*ones(1,2*l);
        %Gradient of the weighting function
        gD =repmat([ -1/2+1/2*s;0],1,2*l);
         %The gradient of the phase functions
        gmu = repmat([((1/2-1/2*t)*a1+(1/2+1/2*t)*a2)/(1-s)+((1/2-1/2*t)*a1+(1/2+1/2*t)*a2)*(1+s)/(1-s)^2;...
                                                                      (-1/2*a1+1/2*a2)*(1+s)/(1-s)],1,2*l);    
        %The Jacobian of the transformation
        J = [ (-x1-x2+x3+x4+x1*t-x2*t-x3*t+x4*t)/(s-1)^2, -1/2*(2*x1*s-2*x2*s-x3-x3*s+x4+x4*s)/(s-1);...
              (-y1-y2+y3+y4+y1*t-y2*t-y3*t+y4*t)/(s-1)^2, -1/2*(2*y1*s-2*y2*s-y3-y3*s+y4+y4*s)/(s-1)];
        %Inverse transpose of the Jacobian
        JTi2 = inv(J.'*J);
        %Determinant of the Jacobian matrix
        detJ = abs(det(J));
        %Calculate matrices
        Ke = Ke + ws(c1)*wt(c2)*((gN.'*JTi2*gD*diag(N)+gN.'*JTi2*gN*diag(D)).'*detJ);
        Ce = Ce + ws(c1)*wt(c2)*((dot(gN.'*JTi2,gmu.',2)*(D.*N)-N.'*N.*(gmu.'*JTi2*gD)-N.'*D.*(gmu.'*JTi2*gN))*detJ).'; 
        Me = Me + ws(c1)*wt(c2)*((N.'*(D.*N)-(diag(N)*dot(gmu.'*JTi2,gmu.',2)*(D.*N)))*detJ).';
    end
end

% For testing it can be useful to calculate by elements
% Calculate by elements
%         for c3=1:2*l
%             for c4 = 1:2*l
%                 Kee(c3,c4) = Kee(c3,c4) + w(c1)*w(c2)*( D(c3)*gN(:,c3).'*JTi2*gN(:,c4)+...
%                                                       N(c3)*gD(:,c3).' *JTi2*gN(:,c4))*detJ;
%                 Cee(c3,c4) = Cee(c3,c4) + w(c1)*w(c2)*( gmu(:,c3).'*JTi2*gN(:,c4)*N(c3)*D(c3)-...
%                                                         gN(:,c3).'*JTi2*gmu(:,c4)*D(c3)*N(c4)-...
%                                                         N(c3)*gD(:,c3).'*JTi2*gmu(:,c4)*N(c4))*detJ;
%                 Mee(c3,c4) = Mee(c3,c4) + w(c1)*w(c2)*( N(c3)*D(c3)*N(c4)-N(c3)*D(c3)*N(c4)*gmu(:,c3).'*JTi2*gmu(:,c4))*detJ;
%             end
%         end
