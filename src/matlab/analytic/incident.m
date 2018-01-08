function [p, q] = incident(type, varargin)
%INCIDENT Incident pressure and velocity wave field
%   [p, q] = incident('point', r0, r, n, k, symm)
%   [p, q] = incident('line', r0, r, n, k, symm)
%   [p, q] = incident('plane', dir, r, n, k, symm)

% Peter Fiala
% 2009-2013

r = varargin{2};
n = varargin{3};
k = varargin{4};
N = size(r,1);

switch lower(type)
    case 'point'
        r0 = varargin{1};
        M = size(r0,1);
        if length(varargin) > 4
            symm = varargin{5};
        else
            symm = 0;
        end
        p = zeros(N,M);
        if nargout == 2
            q = zeros(N,M);
        end
        for m = 1 : size(r0,1)
            dvec = r - repmat(r0(m,:), N, 1);
            d = sqrt(dot(dvec, dvec, 2));
            p(:,m) = p(:,m) + exp(-1i*k*d)./ d / (4*pi);
            if nargout == 2
                rdn = dot(dvec, n, 2) ./ d;
                q(:,m) = q(:,m) -(1+1i*k*d)./d .* p(:,m) .* rdn;
            end
        end
        if symm
            [p2, q2] = incident('point', r0*diag([1 1 -1]), r, n, k);
            p = p + symm*p2;
            q = q + symm*q2;
        end
    case 'line'
        r0 = varargin{1};
        M = size(r0,1);

        p = zeros(N,M);
        if nargout == 2 
            q = zeros(N,M);
        end
        for m = 1 : size(r0,1)
            dvec = r(:,1:2) - repmat(r0(m,1:2), N, 1);
            d = sqrt(dot(dvec, dvec, 2));
            p(:,m) = p(:,m) -1i/4 * besselh(0,2,k*d);
            if nargout == 2
                rdn = dot(dvec, n(:,1:2), 2) ./ d;
                q(:,m) = q(:,m) + 1i*k/4 * besselh(1,2,k*d) .* rdn;
            end
        end
    case 'plane'
        dir = varargin{1};
        dir = dir/norm(dir);
        if length(varargin) > 4
            symm = varargin{5};
        else
            symm = 0;
        end
        kk = dir * k;
        p = exp(-1i*(r*kk(:)));
        if nargout == 2
            q = -1i*(n*kk(:)) .* p;
        end
        if symm
            [p2, q2] = incident('plane', dir*diag([1 1 -1]), r, n, k);
            p = p + symm*p2;
            q = q + symm*q2;
        end
    otherwise
end % of switch

end % of function
