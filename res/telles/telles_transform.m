function [xi, w] = telles_transform(eta, wi, eta_bar, d)
%TELLES_TRANSFORM Third order polynomial quadrature transform by Telles
%   The Telles transformation transform a quadrature into a new one such
%   that certain singularities can be cancelled or near-singularities can
%   be integrated more accurately using numerical integration.
%
%   [XI, W] = TELLES_TRANSFORM(ETA, WI, ETA_BAR, D) Transforms the original
%       quadrature defined by the base points ETA and weights WI to a new
%       quadrature XI, W. ETA should be of type N x DIM and WI must be of 
%       type N x 1. XI and W will have the same size and ETA and WI, 
%       respectively. 
%       The transformation concentrates the base points around the singular 
%       point defined by ETA_BAR. Note that ETA_BAR is not necessarily 
%       inside the domain of integration. The optional parameter D is the 
%       distance from the singularity in the other direction, also in the 
%       coordinate system of ETA. If no D distance is given, D = 0 is
%       assumed.
%
%   References:
%   [1] J.C.F. Telles: A self-adaptive co-ordinate transformation for
%       efficient numerical evaluation of general boundary element
%       integrals. International Journal for numerical methods in
%       engineering. Vol. 24, pp. 959-973 (1987)
%       DOI: 10.1002/nme.1620240509

% Distance parameter
if nargin < 4
    d = 0;
end

% Distance dependent parametrization
if d < 0.05
    r_bar = 0;
elseif d <= 1.3
    r_bar = 0.85 + 0.24 * log(d);
elseif d <= 3.618
    r_bar = 0.893 + 0.0832 * log(d);
else
    r_bar = 1;
end

n_dim = size(eta, 2);
xi = zeros(size(eta));
w = wi;

% Transform along each dimension
for i_dim = 1 : n_dim
    e = eta_bar(i_dim);         % Shortcut notation 
    f = (1 + 2*r_bar);          % Helper quantity
    % As defined in eq. (21)
    q = 1/(2*f)*((e*(3-2*r_bar) - 2*e^3/f)*1/f - e);
    p = 1/(3*f^2)*(4*r_bar*(1-r_bar) + 3*(1 - e^2));
    
    s = sqrt(q^2 + p^3);
    % The singular point in gamma
    gamma_bar = nthroot(-q+s, 3) + nthroot(-q-s, 3) + e/f;
    
    Q = 1 + 3*gamma_bar^2;
    % Polynomial coefficients as defined in eq. (21)
    a = (1 - r_bar) / Q;
    b = -3*(1 - r_bar) * gamma_bar / Q;
    c = (3*gamma_bar^2 + r_bar) / Q;
    d = -b;

    % Evaluate polynomial transform
    xi(:, i_dim) = a * eta(:, i_dim).^3 + b * eta(:, i_dim).^2 + ...
                   c * eta(:, i_dim) + d; 
    
    % Jacobian of the polynomial transform           
    jac = 3*a*eta(:, i_dim).^2 + 2*b*eta(:, i_dim) + c;
    w = w .* jac;
end

end % of function telles_transform

