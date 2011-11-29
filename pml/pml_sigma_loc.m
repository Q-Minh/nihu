% Absorbing function
function [s] = pml_sigma_loc(x, c, L)

% x - coordinates
% c - sound velocity
% L - PML length

s = c./(1-x)/L;
end