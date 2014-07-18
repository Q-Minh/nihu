function [Ke,Ce,Me] = iemat34m(M,z,ngauss)
%IEMAT34M Computes element matrices for 3D hexahedric infinite element.
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
% See also: iemap3D4, lagpolys, gaussquad, iepproc34m

% Peter Rucz
% 2009 november

%Create lagrangian shape functions 
l = length(z);
L = lagpolys(z);
dL = cell(l,1);
% Scale the shape functions
for c1 = 1:length(z);
    L{c1} = L{c1}/polyval(L{c1},z(c1))/(1-z(c1));
    temp = L{c1}.*(length(L{c1})-1:-1:0);
    dL{c1} = temp(1:end-1);
end

% Coordinates of the mapping
x1=M(1,1); y1=M(1,2); z1=M(1,3);
x2=M(2,1); y2=M(2,2); z2=M(2,3);
x3=M(3,1); y3=M(3,2); z3=M(3,3);
x4=M(4,1); y4=M(4,2); z4=M(4,3);
x5=M(5,1); y5=M(5,2); z5=M(5,3);
x6=M(6,1); y6=M(6,2); z6=M(6,3);
x7=M(7,1); y7=M(7,2); z7=M(7,3);
x8=M(8,1); y8=M(8,2); z8=M(8,3);

%Calculate distance of nodes from location of virtual source points
a1 = norm(M(1,:)-M(5,:));
a2 = norm(M(2,:)-M(6,:));
a3 = norm(M(3,:)-M(7,:));
a4 = norm(M(4,:)-M(8,:));

%Gaussian quadrature base points and weights
[ps,ws] = gaussquad(ngauss(1));
[pt,wt] = gaussquad(ngauss(2));
[pu,wu] = gaussquad(ngauss(3));

Ke = zeros(4*l);
Ce = zeros(4*l);
Me = zeros(4*l);

N = zeros(1,4*l);
gN = zeros(3,4*l);

for c1 = 1:ngauss(1)
    for c2 = 1:ngauss(2)
        for c3 = 1:ngauss(3)
            s = ps(c1);
            t = pt(c2);
            u = pu(c3);
            %Calculate shape functions and their gradients
            for c4 = 1:l
                % Shape functions
                N(4*(c4-1)+1) = 1/4*(1-s)*polyval(L{c4},s)*(1-t)*(1-u);
                N(4*(c4-1)+2) = 1/4*(1-s)*polyval(L{c4},s)*(1+t)*(1-u);
                N(4*(c4-1)+3) = 1/4*(1-s)*polyval(L{c4},s)*(1+t)*(1+u);
                N(4*(c4-1)+4) = 1/4*(1-s)*polyval(L{c4},s)*(1-t)*(1+u);
                % Calculate gradients
                % d/ds
                gN(1,4*(c4-1)+1) = 1/4*(1-t)*(1-u)*((-1)*polyval(L{c4},s)+(1-s)*polyval(dL{c4},s));
                gN(1,4*(c4-1)+2) = 1/4*(1+t)*(1-u)*((-1)*polyval(L{c4},s)+(1-s)*polyval(dL{c4},s));
                gN(1,4*(c4-1)+3) = 1/4*(1+t)*(1+u)*((-1)*polyval(L{c4},s)+(1-s)*polyval(dL{c4},s));
                gN(1,4*(c4-1)+4) = 1/4*(1-t)*(1+u)*((-1)*polyval(L{c4},s)+(1-s)*polyval(dL{c4},s));
                % d/dt
                gN(2,4*(c4-1)+1) = 1/4*(1-s)*polyval(L{c4},s)*(1-u)*(-1);
                gN(2,4*(c4-1)+2) = 1/4*(1-s)*polyval(L{c4},s)*(1-u)*( 1);
                gN(2,4*(c4-1)+3) = 1/4*(1-s)*polyval(L{c4},s)*(1+u)*( 1);
                gN(2,4*(c4-1)+4) = 1/4*(1-s)*polyval(L{c4},s)*(1+u)*(-1);
                % d/du
                gN(3,4*(c4-1)+1) = 1/4*(1-s)*polyval(L{c4},s)*(1-t)*(-1);
                gN(3,4*(c4-1)+2) = 1/4*(1-s)*polyval(L{c4},s)*(1+t)*(-1);
                gN(3,4*(c4-1)+3) = 1/4*(1-s)*polyval(L{c4},s)*(1+t)*( 1);
                gN(3,4*(c4-1)+4) = 1/4*(1-s)*polyval(L{c4},s)*(1-t)*( 1);

            end
            %Weighting function
            D = (1-s)^2/4*ones(1,4*l);
            %Gradient of the weighting function
            gD =repmat([ -1/2+1/2*s;0;0],1,4*l);
            %The gradient of the phase functions
            gmu = repmat([ 1/2*(a1-a1*u-a1*t+a1*u*t+a2-a2*u+a2*t-a2*u*t+a3+a3*u+a3*t+a3*u*t+a4+a4*u-a4*t-a4*u*t)/(s-1)^2;...
                          -1/4*(-a1+a1*u+a2-a2*u+a3+a3*u-a4-a4*u)*(1+s)/(s-1);...
                          -1/4*(-a1+a1*t-a2-a2*t+a3+a3*t+a4-a4*t)*(1+s)/(s-1)],1,4*l);
            %The Jacobian of the transformation
            J =  [  1/2*(x5-x4+x6+x7+x8-x1*u*t-x3*u*t-x8*u*t+x5*u*t-x1-x2-x3+x2*u*t+x4*u*t+x1*t+x1*u-x2*t+x2*u-x3*t-x3*u+x4*t-x4*u-x5*t-x5*u+x6*t-x6*u+x7*t+x7*u-x8*t+x8*u-x6*u*t+x7*u*t)/(s-1)^2,   -1/4*(2*x1*s-2*x1*u*s-2*x2*s+2*x2*u*s-2*x3*s-2*x3*u*s+2*x4*s+2*x4*u*s-x5-x5*s+x5*u+x5*u*s+x6+x6*s-x6*u-x6*u*s+x7+x7*s+x7*u+x7*u*s-x8-x8*s-x8*u-x8*u*s)/(s-1),  -1/4*(2*x1*s-2*x1*t*s+2*x2*s+2*x2*t*s-2*x3*s-2*x3*t*s-2*x4*s+2*x4*t*s-x5-x5*s+x5*t+x5*t*s-x6-x6*s-x6*t-x6*t*s+x7+x7*s+x7*t+x7*t*s+x8+x8*s-x8*t-x8*t*s)/(s-1);...
                1/2*(-y4+y5+y6+y7+y8+y5*u*t+y4*u*t-y6*u*t+y2*u*t-y3*u*t-y1-y2-y3+y1*t+y1*u-y2*t+y2*u-y3*t-y3*u+y4*t-y4*u-y5*t-y5*u-y1*u*t+y7*u*t-y8*u*t+y6*t-y6*u+y7*t+y7*u-y8*t+y8*u)/(s-1)^2,  -1/4*(2*y1*s-2*y1*u*s-2*y2*s+2*y2*u*s-2*y3*s-2*y3*u*s+2*y4*s+2*y4*u*s-y5-y5*s+y5*u+y5*u*s+y6+y6*s-y6*u-y6*u*s+y7+y7*s+y7*u+y7*u*s-y8-y8*s-y8*u-y8*u*s)/(s-1),  -1/4*(2*y1*s-2*y1*t*s+2*y2*s+2*y2*t*s-2*y3*s-2*y3*t*s-2*y4*s+2*y4*t*s-y5-y5*s+y5*t+y5*t*s-y6-y6*s-y6*t-y6*t*s+y7+y7*s+y7*t+y7*t*s+y8+y8*s-y8*t-y8*t*s)/(s-1);...
                1/2*(-z4+z5+z6+z7+z8-z1*u*t+z5*u*t+z4*u*t-z3*u*t-z8*u*t+z2*u*t-z1-z2-z3-z6*u*t+z7*u*t+z1*t+z1*u-z2*t+z2*u-z3*t-z3*u+z4*t-z4*u-z5*t-z5*u+z6*t-z6*u+z7*t+z7*u-z8*t+z8*u)/(s-1)^2,  -1/4*(2*z1*s-2*z1*u*s-2*z2*s+2*z2*u*s-2*z3*s-2*z3*u*s+2*z4*s+2*z4*u*s-z5-z5*s+z5*u+z5*u*s+z6+z6*s-z6*u-z6*u*s+z7+z7*s+z7*u+z7*u*s-z8-z8*s-z8*u-z8*u*s)/(s-1),  -1/4*(2*z1*s-2*z1*t*s+2*z2*s+2*z2*t*s-2*z3*s-2*z3*t*s-2*z4*s+2*z4*t*s-z5-z5*s+z5*t+z5*t*s-z6-z6*s-z6*t-z6*t*s+z7+z7*s+z7*t+z7*t*s+z8+z8*s-z8*t-z8*t*s)/(s-1)];
            %Inverse transpose of the Jacobian
            JTi2 = inv(J.'*J);
            %Determinant of the Jacobian matrix
            detJ = abs(det(J));
            %Calculate matrices
            Ke = Ke + ws(c1)*wt(c2)*wu(c3)*((gN.'*JTi2*gD*diag(N)+gN.'*JTi2*gN*diag(D)).'*detJ);
            Ce = Ce + ws(c1)*wt(c2)*wu(c3)*((dot(gN.'*JTi2,gmu.',2)*(D.*N)-N.'*N.*(gmu.'*JTi2*gD)-N.'*D.*(gmu.'*JTi2*gN))*detJ).';
            Me = Me + ws(c1)*wt(c2)*wu(c3)*((N.'*(D.*N)-(diag(N)*dot(gmu.'*JTi2,gmu.',2)*(D.*N)))*detJ).';
        end
    end
end
