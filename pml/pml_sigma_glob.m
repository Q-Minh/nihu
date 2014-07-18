% Absorbing function for gamma_x,y,z
function [s] = pml_sigma_glob(x, a, c)

% x - coordinates
% c - sound velocity
% L - PML length

s = c./abs(a-x);