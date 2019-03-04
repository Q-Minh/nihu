function [xi, w] = telles_transform(eta, wi, eta_bar, d, alpha)
%TELLES_TRANSFORM Third order polynomial quadrature transform by Telles
%   The Telles transformation transform a quadrature into a new one such
%   that certain singularities can be cancelled or near-singularities can
%   be integrated more accurately using numerical integration.
%
%   [XI, W] = TELLES_TRANSFORM(ETA, WI, ETA_BAR, D, ALPHA) Transforms the 
%       original quadrature defined by the base points ETA and weights WI 
%       to a new quadrature XI, W. ETA should be of type N x DIM and WI 
%       must be of type N x 1. XI and W will have the same size and ETA and 
%       WI, respectively. 
%       The transformation concentrates the base points around the singular 
%       point defined by ETA_BAR. Note that ETA_BAR is not necessarily 
%       inside the domain of integration. The optional parameter D is the 
%       distance from the singularity in the other direction, also in the 
%       coordinate system of ETA. If no D distance is given, D = 0 is
%       assumed.
%       The optional parameter ALPHA optimizes the choice of the reference
%       distance for the transformation based on the singularity type:
%       O(1/r^ALPHA), (ALPHA = 0 means singularity of log(r) type). If
%       ALPHA is not set, ALPHA = 0 assumed.
%
%   References:
%   [1] J.C.F. Telles: A self-adaptive co-ordinate transformation for
%       efficient numerical evaluation of general boundary element
%       integrals. International Journal for numerical methods in
%       engineering. Vol. 24, pp. 959-973 (1987)
%       DOI: 10.1002/nme.1620240509
%
%   [2] J.C.F. Telles, R.F. Oliveira: Third degree polynomial
%       transformation for boundary element integrals: Further
%       improvements. Engineering Analysis with boundary elements. Vol. 13,
%       pp. 135-141 (1994)
%       DOI: 10.1016/0955-7997(94)90016-7
%      

% Tabulated data from [2], Table 1.
r_bar_ref = [
 % D        log(r)      1/r         1/r^2      1/r^3      
 0.00       0           0           0          0
 0.01       0.06915407  0           0          0 
 0.016667   0.09831843  0.05747304  0          0
 0.05       0.4483645   0.2050263   0.1377359  0.1118472
 0.10       0.5287795   0.3242594   0.2239878  0.1848561
 0.20       0.6446450   0.4830719   0.3731488  0.3028686 
 0.30       0.7215894   0.5895026   0.4916305  0.4179902
 0.50       0.8181035   0.7270662   0.6519612  0.5899927
 0.70       0.8747554   0.8094936   0.7523058  0.7024217
 0.90       0.9101131   0.8621070   0.8181204  0.7783351
 1.50       0.9588171   0.9374435   0.9162750  0.8956607
 2.50       0.9821442   0.9735865   0.9652033  0.9569454
 3.50       0.9899879   0.9852044   0.9805508  0.9760106
 4.50       0.9936338   0.9905608   0.9875566  0.9846170
14.00       0.9992856   0.9989301   0.9985758  0.9982227
15.00       1.0         1.0         1.0        1.0 
16.00       1.0         1.0         1.0        1.0 
];

if nargin < 4
    d = 0;
end


if nargin < 5
    alpha = 0;
end

r_bar = interp1(r_bar_ref(:,1), r_bar_ref(:, alpha+2), d, 'linear', 'extrap');

% Distance parameter

% % Distance dependent parametrization
% if d < 0.05
%     r_bar = 0;
% elseif d <= 1.3
%     r_bar = 0.85 + 0.24 * log(d);
% elseif d <= 3.618
%     r_bar = 0.893 + 0.0832 * log(d);
% else
%     r_bar = 1;
% end

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

