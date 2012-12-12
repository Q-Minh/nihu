function [vtil, kz] = trans_inf_invar(vz, z)
%TRANS_INF_INVAR Compute wavenumber domain excitation of the FE model
%   [VTIL KZ] = TRANS_INF_INVAR(VZ, Z) Computes the wavenumber domain
%   excitation vtil from the spatial domain excitation vz and the spatial
%   cooridnates z

% Peter Fiala
% 2008 June

nz = size(vz,2)-1;
Z = max(z)-min(z);
vtil = Z*fftshift(ifft(vz(:,1:end-1), [], 2),2);
kz = linspace(-nz*pi/Z, nz*pi/Z, nz+1);
kz = kz(1:end-1);