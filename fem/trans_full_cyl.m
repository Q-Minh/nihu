function [vtil, m] = trans_full_cyl(vphi, nPhi)
%TRANS_FULL_CYL Compute wavenumber domain excitation of the FE model
%   [VTIL M] = TRANS_FULL_CYL(VPHI, NPHI) Computes the wavenumber domain
%   excitation vtil from the spatial domain excitation vphi and the angular
%   cooridnates phi

% Peter Fiala
% 2008 June

vtil = fftshift(fft(vphi, [], 2), 2)/nPhi;
m = -nPhi/2 : nPhi/2-1;