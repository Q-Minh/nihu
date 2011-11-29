function [xi, err, it] = invmap134(X,M,tol,ctol)
% INVMAP38 Inverse numerical isoparametric mapping for hexa elements.
%   [XI, ERR, IT] = INVMAP38(X, M, TOL, CTOL) returns XI = [xi, eta, zeta]
%       local coordinates for the point with global coordinates
%       M = (x, y, z) in the hexa element given by the node coordinates X
%       in a 8x3 matrix. TOL is the relative tolerance for the location
%       error. CTOL is the tolerance parameter for the convergence in the
%       iteration.

% Peter Rucz
% Last modified: 2010.04.07.

%% Argument check and default parameters
%error(nargchk(2,4,nargin,'struct'));
%error(nargoutchk(1,3,nargout,'struct'));

% Check if inside
T1 = sum(X(1:4,:),1)/4;
T2 = sum(X(5:8,:),1)/4;
T1M = M-T1;
N12 = norm(T1-T2);
NT = (T2-T1)/N12;
% dev1 = max(sqrt(dot(repmat(T1,4,1)-X(1:4,:),repmat(T1,4,1)-X(1:4,:),2)));
% dev2 = max(sqrt(dot(repmat(T2,4,1)-X(5:8,:),repmat(T2,4,1)-X(5:8,:),2)));
dev1 = max([T1;T1;T1;T1]-X(1:4,:),[],1);
dev2 = max([T2;T2;T2;T2]-X(5:8,:),[],1);

l = T1M*NT.';
DM = norm(M-(T1+l*NT));
if(DM > dev1 + (dev2 - dev1) * l/N12)
    %disp('Bambo!');
    xi = [inf inf inf];
    return
end

if nargin < 4
    ctol = 1e2;
end
if nargin < 3
    tol = 1e-5;
end
coefftol = 1e-6;
maxit = 100;
Dtol = tol*norm(X(1,:)-X(2,:));

% Coefficients for shape functions
Ncoeff = 1/2*[
    1     0    -1     0     -1    0    1     0
   -1     0    -1     0      1    0    1     0
    1     0     1     0      1    0    1     0
   -1     0     1     0     -1    0    1     0;
    -1/2*[
    1     1    -1    -1    -1    -1     1     1
   -1    -1    -1    -1     1     1     1     1
    1     1     1     1     1     1     1     1
   -1    -1     1     1    -1    -1     1     1]
    ];

cconst = Ncoeff(:,[2 6 4 8]);
czeta  = Ncoeff(:,[1 5 3 7]);
cc = [cconst, czeta];
converged = false;



% V1 = X - repmat(X(1,:),8,1);
% V7 = X - repmat(X(7,:),8,1);
% C1 = cross(V1([2 2 3 4 4 5 2 2 5],:),V1([4 3 4 5 8 8 5 6 6],:));
% C7 = cross(V7([3 2 2 3 3 4],:),V7([6 3 6 8 4 8],:));
% N1 = C1 .* repmat(sign(dot(V1([5 5 5 2 2 2 4 4 4],:),C1,2)),1,3);
% N7 = C7 .* repmat(sign(dot(V7([8 8 8 6 6 6],:),C7,2)),1,3);
% d1   = M - X(1,:);
% d7   = M - X(7,:);
% Dt = reshape(dot([N1 ; N7],[repmat(d1,9,1); repmat(d7,6,1)],2),3,5).';
%
% if(sum(all(Dt < 0,2)) > 0)
%     %disp('Rambo!');
%     xi = [inf inf inf];
%     return
% end

%% Setup initial iteration
for p = 1:8
    d = M - X(p,:);
    n1 = d*[0 0 0; 0 0 1; 0 -1 0];
    n2 = d*[0 1 0; -1 0 0; 0 0 0];
    if sum(n1~=0) < 2 || sum(n2~=0) < 2
        % go to next vertex
        continue;
    end
  
    ccc = -cc.'*(X*[n1.' n2.']);
    % Modified here + -> -
    ccc(4,:) = ccc(4,:) - X(p,:)*[n1.' n2.'];
    ccc(8,:) = ccc(8,:) + X(p,:)*[n1.' n2.'];
        
    det17 = det(ccc([1 7],:));
    det35 = det(ccc([3 5],:));
    det57 = det(ccc([5 7],:));
    det18 = det(ccc([1 8],:));
    det27 = det(ccc([2 7],:));
    det36 = det(ccc([3 6],:));
    det45 = det(ccc([4 5],:));
    det58 = det(ccc([5 8],:));
    det67 = det(ccc([6 7],:));
    det28 = det(ccc([2 8],:));
    det46 = det(ccc([4 6],:));
    det68 = det(ccc([6 8],:));
    
   % Guess initial zeta
    zeta = 0;
    dzeta = 0;
    D = [inf inf inf];  % distance between guess location and M
    it = 0;
    %% The main iteration
    nextvertex = false;
    while(norm(D) > Dtol)
        it = it + 1;
        if it > maxit
            nextvertex = true;
            break;
        end
        zeta = zeta+dzeta;
        
        ccs = ccc(1:4,:) + ccc(5:8,:)*zeta;
        % coefficients for quadratic solution
        a = det(ccc([1 3],:) + ccc([5 7],:)*zeta);
        b = det(ccc([1 4],:) + ccc([5 8],:)*zeta)+...
            det(ccc([2 3],:) + ccc([6 7],:)*zeta);
        c = det(ccc([2 4],:) + ccc([6 8],:)*zeta);
        
        md = max(abs([a b c]));
        
        if (abs(a) / md < coefftol && abs(b)/md < coefftol) 
            nextvertex = true;
            break
        end
        
        % solve a + bx1 + cx2 + dx1x2 = 0
        xi = zeros(2);
        xi(:,1) = quadratic(a, b, c);
        xi(:,2) = -((ccc(4,2) + ccc(8,2)*zeta) + (ccc(3,2)+ccc(7,2)*zeta)*xi(:,1)) ./...
               ((ccc(2,2)+ccc(6,2)*zeta) + (ccc(1,2)+ccc(5,2)*zeta)*xi(:,1));
        xi = xi(all(abs(xi) < 1e6,2),:);
        if(size(xi,1)>1)
            if(norm(xi(1,:)) < norm(xi(2,:)))
                xi = xi(1,:);
                sel = 0;
            else
                xi = xi(2,:);
                sel = 1;
            end
        else
            sel = 0;
        end
        eta = xi(2);
        xi = xi(1);
      
        % da / dzeta
        dadzeta = det17 - det35 + 2*zeta*det57;
        dbdzeta = det18 + det27 - det36 - det45 + 2*zeta*(det58+det67);
        dcdzeta = det28 - det46 + 2*zeta*det68;
        
        % dxi / dzeta
        if (a ~= 0)
            dxidzeta = (b/2/a^2 + (-1)^sel * ((2*a*c - b^2)/(2*a^2*(b^2 - 4*a*c)^(1/2))))*dadzeta + ...
                (-1/2/a  + (-1)^sel * (   b/(sqrt(b^2-4*a*c)))/2/a)*dbdzeta + ...
                + (-1)^sel * (-2*a/(sqrt(b^2-4*a*c))/2/a)*dcdzeta;
        else
            dxidzeta = -1/b*dcdzeta + c/b^2*dbdzeta;
        end
        % deta / dxi
        detadxi = (-ccs(3,2)*(ccs(1,2)*xi+...
            ccs(2,2)+ccs(1,2)*(ccs(3,2)*xi+ccs(4,2))))...
            /(ccs(1,2)*xi+ccs(2,2))^2;
        if(isnan(detadxi))
            detadxi = 1;
        end
        % deta / dzeta
        detadzeta = detadxi * dxidzeta;
        
        % Calculate derivatives
        Nxi = [xi*eta*zeta xi*eta xi*zeta xi eta*zeta eta zeta 1]*Ncoeff.'/(zeta-1);
        dNxi = [[eta*zeta eta zeta  1   0  0 0 0;
                 xi*zeta  xi    0  0 zeta 1 0 0]/(zeta-1);
                -[xi*eta   xi*eta xi xi eta eta 1 1]/(zeta-1)^2]*Ncoeff.';
        dxdzetac = [dxidzeta detadzeta 1] * dNxi * X;
        
        % The guessed location
        Mhat = Nxi*X;
        Dprev = D;
        D = M - Mhat;
        
        if(ctol*norm(Dprev) < norm(D))
            nextvertex = true;
            break
        else
            dzeta = D*dxdzetac.'/(dxdzetac*dxdzetac.');
            if (zeta + dzeta >= 1)
                dzeta = (1 - zeta)/ 2;
            end
        end
    end
    if ~nextvertex
        converged = true;
        break;
    end
end


if (~converged)
    xi = [inf inf inf];
    disp('Not converged');
else
   % disp(['Converged, iterations = ',num2str(it)]);
    xi = real([zeta xi eta]);
    err = real(D);
end
%disp(p);
end
