function [xi, w] = telles_tria(n, eta_bar, d)
%TELLES_TRIA Summary of this function goes here
%   Detailed explanation goes here

% Create a quadrature over a line element as a base
[xi0, w0] = gaussquad(n);

% Re-transform the singular point in xi_1 from [0, 1] to [-1, 1]
xi1_bar = 2*eta_bar(1) -1;
% The jacobian of the re-transformation in xi_1 
jac1 = 2;
[xi_t, w_xi_t] = telles_transform(xi0, w0, xi1_bar, d * jac1);

% Preallocate for base points and weights
xi = zeros(n*n, 2);
w  = zeros(n*n, 1);

for i = 1 : n
    idx = (i-1)*n + (1:n);
    % Actual coordinate in xi1
    xi1 = (1 + xi_t(i))/2;
    % Re-transform the singular point in xi_2 from [0, 1-xi1] to [-1, 1]
    xi2_bar = 2*eta_bar(2) / (1 - xi1) - 1;
    % The jacobian of the re-transformation in xi_2
    jac2 = 2 / (1 - xi1);
    % Transform the line quadrature 
    [eta_t, w_eta_t] = telles_transform(xi0, w0, xi2_bar, d*jac2);
    % Fill quadrature base points using scaling
    xi(idx, :) = [ xi1 * ones(n, 1), (1-xi1)*(eta_t+1)/2 ];
    % Fill weights using the jacobians as scaling
    w(idx) = (w_xi_t(i) / jac1) * (w_eta_t / jac2);
end

end

