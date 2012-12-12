function [phat] = itrans_full_cyl(ptil, m)
%ITRANS_FULL_CYL  Inverse modal expansion of response
%   PHAT = ITRANS_FULL_CYL(PTIL, M) Transforms the modal coordinates PTIL
%   of the response into the spatial domain PHAT.
%   PTIL : modal coordinates of response (nDof x nM x nFreq matrix)
%   M    : nM element vector
%   PHAT : nDOF x (nM+1) x nFreq matrix
%
% See also: trans_full_cyl, solve_full_cylindrical_direct,
% solve_full_cylindrical_modal

% Peter Fiala
% 2008 June

phat = length(m) * ifft(fftshift(ptil,2),[],2);