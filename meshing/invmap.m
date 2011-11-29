function xi = invmap(xnod, x, type)
%INVMAP Inverse isoparametric mapping
%  XI = INVMAP(XNOD, X, TYPE) returns the isoparametric mapping space
%  coordinates XI of a location X in an element defined by the nodal
%  coordinates XNOD.

%   Copyright 2008-2010 P. Fiala
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

%% Argument check
error(nargchk(3,4,nargin,'struct'));

%% Inverse mapping
switch type
    case 12 % line element
        xi = (2*x - (xnod(1)+xnod(2))) / (xnod(2)-xnod(1));
    case 23 % TRIA element
        xnod = xnod(:,1:2);
        D = [
            -1 1 0;
            -1 0 1
            ];
        xi = (x(1:2) - xnod(1,:)) / (D * xnod);
    case 24 % QUAD element
        a = -x + [1 1 1 1]/4 * xnod;
        b = [-1 1 1 -1]/4 * xnod;
        c = [-1 -1 1 1]/4 * xnod;
        d = [1 -1 1 -1]/4 * xnod;
        if norm(d)/norm(b) < 1e-6 % parallelogram element
            xi = -a / [b;c];
        else % general quad element
            % solve a + bx1 + cx2 + dx1x2 = 0
            xi(:,1) = quadratic(det([d; b]), det([c; b])+det([d; a]),...
                det([c; a]));
            xi(:,2) = -(a(2) + b(2)*xi(:,1)) ./ (c(2) + d(2)*xi(:,1));
            xi = xi(all(abs(xi) < 1e3,2),:);
        end
    case 34 % TETRA element
        D = [
            -1 1 0 0;
            -1 0 1 0
            -1 0 0 1
            ];
        xi = (x - xnod(1,:)) / (D * xnod);
    case 38 % BRICK element
        xi = invmap38b(xnod, x);
    case 122 %2D quad infinite element
        a =  2*x + [0 0 -1 -1] * xnod;
        b = -2*x + [2 2 -1 -1] * xnod;
        c = [0 0  1 -1] * xnod;
        d = [-2 2 1 -1] * xnod;
        % solve a + bx1 + cx2 + dx1x2 = 0
        if(d == 0)
            xi = (-[b.',c.']\a.').';
        else
            xi(:,1) = quadratic(det([d; b]), det([c; b])+det([d; a]),...
                det([c; a]));
            xi(:,2) = -(a(2) + b(2)*xi(:,1)) ./ (c(2) + d(2)*xi(:,1));
            xi = xi(all(abs(xi) < 1e3,2),:);
        end
    case 133 %3D penta infinite element
        xi = invmap133b(xnod, x);
    case 134 %3D hexa infinite element
        xi = invmap134(xnod, x);
    otherwise
        error('nihu:invmap:TODO',...
            'Inverse mapping not yet implemented for type %d', type);
end
end

function x = quadratic(a, b, c)
%QUADRATIC  Solve quadratic polynomial equation ax^2 + bx + c = 0

if a ~= 0
    det = sqrt(b^2-4*a*c);
    x(1) = (-b + det) / (2*a);
    x(2) = (-b - det) / (2*a);
else
    x(1) = -c/b;
    x(2) = Inf;
end
end