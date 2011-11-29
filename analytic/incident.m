function [p, q] = incident(type, varargin)
%INCIDENT Incident pressure and velocity wave field
%   [p, q] = incident('point', r0, r, n, k, symm)
%   [p, q] = incident('plane', dir, r, n, k, symm)

% Peter Fiala
% 2009

switch lower(type)
    case 'point'
        r0 = varargin{1};
        r = varargin{2};
        n = varargin{3};
        k = varargin{4};
        if length(varargin) > 4
            symm = varargin{5};
        else
            symm = 0;
        end
        N = size(r,1);
        p = zeros(N,1);
        q = zeros(N,1);
        for m = 1 : size(r0,1)
            dvec = r - repmat(r0(m,:), N, 1);
            d = sqrt(dot(dvec, dvec, 2));
            dvec = dvec ./ repmat(d, 1, 3);
            p = p + exp(-1i*k*d)./ d;
            if nargout == 2 && ~isempty(n)
                q = q -(1+1i*k*d)./d .* p .* dot(dvec, n, 2);
            end
        end
        if symm
            [p2, q2] = incident('point', r0*diag([1 1 -1]), r, n, k);
            p = p + symm*p2;
            q = q + symm*q2;
        end
    case 'line'
    case 'plane'
        dir = varargin{1};
        r = varargin{2};
        n = varargin{3};
        k = varargin{4};
        if length(varargin) > 4
            symm = varargin{5};
        else
            symm = 0;
        end
        dir = dir/norm(dir);
        k = dir * k;
        p = exp(-1i*(r*k(:)));
        if nargout == 2
            q = -1i*(n*k(:)) .* p;
        end
        if symm
            [p2, q2] = incident('plane', dir*diag([1 1 -1]), r, n, k);
            p = p + symm*p2;
            q = q + symm*q2;
        end
    otherwise
end
end
