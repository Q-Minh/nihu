function [xi, w] = telles_tria(n, eta_bar, d)
%TELLES_TRIA Summary of this function goes here
%   Detailed explanation goes here

[xi0, w0] = gaussquad(n);

xi_bar = 2*eta_bar(1) -1;
[xi_t, w_xi_t] = telles_transform(xi0, w0, xi_bar, d/(0.5));

xi = [];
w = [];

for i = 1 : size(xi_t, 1)
    xi_bar = 2*eta_bar(2) / (1/2 - xi_t(i)/2) - 1;
    [eta_t, w_eta_t] = telles_transform(xi0, w0, xi_bar, 2*d/(1/2 - xi_t(i)/2));
    
    xi = [xi; [
            repmat(0.5*xi_t(i)+0.5, size(eta_t, 1), 1), ...
          (0.5*eta_t + 0.5)*(0.5 - 0.5*xi_t(i))]];
    w = [w; .5 * repmat(w_xi_t(i)*size(eta_t, 1), 1) .* (.5*w_eta_t*(0.5 - 0.5*xi_t(i)))];
end

w = w / n;

end

