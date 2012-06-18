function [p, q] = incident_kx(type, varargin)
%INCIDENT_KX Incident wave field in longitudinal wavenumber domain
%   [p, q] = INCIDENT_KX('point', r0, r, n, k, kx)

switch type
    case 'point'
        r0 = varargin{1};   % source location
        r = varargin{2};    % receiver locations
        n = varargin{3};    % unit normal vectors
        k = varargin{4};    % wave number (om / c)
        kx = varargin{5};   % longitudinal wave number
        
        N = size(r,1);      % number of receivers
        p = zeros(N,1);     % incident pressure
        if nargout > 1
            q = zeros(N,1); % incident velocity
        end
        for m = 1 : size(r0,1)  % for each source
            dvec = r(:,2:3) - repmat(r0(m,2:3), N, 1);   % distance vector
            d = sqrt(dot(dvec, dvec, 2));       % transversal distance
            rdn = dot(dvec, n(:,2:3), 2)./d;    % normal derivative of transversal distance
            kt = -1i*sqrt(kx^2-k^2);            % transversal wave number
            padd = -1i/4 * besselh(0,2, kt*d);
            p = p + padd;   % incident pressure
            if nargin > 1
                q = q - 1i*kx.*n(:,1).*padd +...   % incident velocity
                    1i*kt*rdn/4 .* besselh(1, 2, kt*d);
            end
        end
    otherwise
        error('invalid type parameter: ''%s''', type);
end
end
