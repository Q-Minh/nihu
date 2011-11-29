function [Ke,Ce,Me] = iemat33mj(M,P,ngauss,alpha,beta)
%IEMAT33M Computes element matrices for 3D pentahedric infinite element.
%   M      : Coordinates of the mapping nodes, must be a 6x3 matrix.
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
l = P;

% Coordinates of the mapping
x1=M(1,1); y1=M(1,2); z1=M(1,3);
x2=M(2,1); y2=M(2,2); z2=M(2,3);
x3=M(3,1); y3=M(3,2); z3=M(3,3);
x4=M(4,1); y4=M(4,2); z4=M(4,3);
x5=M(5,1); y5=M(5,2); z5=M(5,3);
x6=M(6,1); y6=M(6,2); z6=M(6,3);

%Calculate distance of nodes from location of virtual source points
a1 = norm(M(1,:)-M(4,:));
a2 = norm(M(2,:)-M(5,:));
a3 = norm(M(3,:)-M(6,:));

%Gaussian quadrature base points and weights
[ps,ws] = gaussquad(ngauss(1));
[ptu,wtu] = dunavant_rule(ngauss(2));
%Correct dunavant rules to match parent domain
ptu(1,:) = -2*ptu(1,:)+1;
ptu(2,:) =  2*ptu(2,:)-1;
wtu      =  2*wtu;

Ke = zeros(3*l);
Ce = zeros(3*l);
Me = zeros(3*l);

N = zeros(1,3*l);
gN = zeros(3,3*l);

for c1 = 1:ngauss(1)
    for c2 = 1:size(ptu,2)
        s = ps(c1);
        t = ptu(1,c2);
        u = ptu(2,c2);
        %Calculate shape functions and their gradients
        for c4 = 1:l
            % Shape functions
            N(3*(c4-1)+1) = 1/4*(1-s)*jacpoly(c4-1,alpha,beta,s)*(1-t);
            N(3*(c4-1)+2) = 1/2*(1-s)*jacpoly(c4-1,alpha,beta,s)*(1-(1+u)/2-(1-t)/2);
            N(3*(c4-1)+3) = 1/4*(1-s)*jacpoly(c4-1,alpha,beta,s)*(1+u);
            % Calculate gradients
            % d/ds
            gN(1,3*(c4-1)+1) =               1/4*(1-t)*((-1)*jacpoly(c4-1,alpha,beta,s)+(1-s)*djacpoly(c4-1,alpha,beta,s,1));
            gN(1,3*(c4-1)+2) = 1/2*(1-(1+u)/2-(1-t)/2)*((-1)*jacpoly(c4-1,alpha,beta,s)+(1-s)*djacpoly(c4-1,alpha,beta,s,1));
            gN(1,3*(c4-1)+3) =               1/4*(1+u)*((-1)*jacpoly(c4-1,alpha,beta,s)+(1-s)*djacpoly(c4-1,alpha,beta,s,1));
            % d/dt
            gN(2,3*(c4-1)+1) = 1/4*(1-s)*jacpoly(c4-1,alpha,beta,s)*(-1);
            gN(2,3*(c4-1)+2) = 1/4*(1-s)*jacpoly(c4-1,alpha,beta,s)*( 1);
            gN(2,3*(c4-1)+3) =                               0;
            % d/du
            gN(3,3*(c4-1)+1) = 0;
            gN(3,3*(c4-1)+2) = 1/4*(1-s)*jacpoly(c4-1,alpha,beta,s)*(-1);
            gN(3,3*(c4-1)+3) = 1/4*(1-s)*jacpoly(c4-1,alpha,beta,s)*( 1);
        end
        %Weighting function
        D = (1-s)^2/4*ones(1,3*l);
        %Gradient of the weighting function
        gD =repmat([ -1/2+1/2*s;0;0],1,3*l);
        %Distance along the inner plate
        % a = (1-t)/2*a1 + (1+u)/2*a2 + (1-(1+u)/2 - (1-t)/2)*a3
        %The gradient of the phase functions
        gmu = repmat([ (a1-a1*t+a2+a2*u-a3*u+a3*t)/(s-1)^2;
            1/2*(a1-a3)*(1+s)/(s-1);
            -1/2*(a2-a3)*(1+s)/(s-1)],1,3*l);
        %The Jacobian of the transformation
        J =[ (-x3*u+x1*t+x2*u-x2*t-x4*t-x5*u+x5*t+x6*u-x1-x3+x4+x6)/(s-1)^2, -1/2*(2*s*x1-2*s*x2-x4-x4*s+x5+x5*s)/(s-1), -1/2*(2*s*x2-2*s*x3-x5-x5*s+x6+x6*s)/(s-1);...
            ( y1*t+y2*u-y3*u-y2*t-y1-y3+y4+y6-y4*t-y5*u+y5*t+y6*u)/(s-1)^2, -1/2*(2*s*y1-2*s*y2-y4-y4*s+y5+y5*s)/(s-1), -1/2*(2*s*y2-2*s*y3-y5-y5*s+y6+y6*s)/(s-1);...
            ( z1*t-z3*u+z2*u-z1-z3+z4+z6-z2*t-z4*t-z5*u+z5*t+z6*u)/(s-1)^2, -1/2*(2*s*z1-2*s*z2-z4-z4*s+z5+z5*s)/(s-1), -1/2*(2*s*z2-2*s*z3-z5-z5*s+z6+z6*s)/(s-1)];
        %Inverse transpose of the Jacobian
        JTi2 = inv(J.'*J);
        %Determinant of the Jacobian matrix
        detJ = abs(det(J));
        %Calculate matrices
        Ke = Ke + ws(c1)*wtu(c2)*((gN.'*JTi2*gD*diag(N)+gN.'*JTi2*gN*diag(D)).'*detJ);
        Ce = Ce + ws(c1)*wtu(c2)*((dot(gN.'*JTi2,gmu.',2)*(D.*N)-N.'*N.*(gmu.'*JTi2*gD)-N.'*D.*(gmu.'*JTi2*gN))*detJ).';
        Me = Me + ws(c1)*wtu(c2)*((N.'*(D.*N)-(diag(N)*dot(gmu.'*JTi2,gmu.',2)*(D.*N)))*detJ).';

    end
end
