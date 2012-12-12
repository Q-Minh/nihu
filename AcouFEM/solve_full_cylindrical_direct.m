function ptil = solve_full_cylindrical_direct(K, M, B, qtil, m, freqvec)
%SOLVE_FULL_CYLINDRICAL_DIRECT Direct solution of cyl symm FE model
%   PTIL = SOLVE_FULL_CYLINDRICAL_DIRECT(K, M, B, QTIL, M, F)
%   solves the cylindrically symmetric FE system in the wavenumber domain.
%   PTIL    : nDOF x nM x nF response in the wavenumber domain
%   K, M, B : nDOF x nDOF FE stiffness and mass matrices
%   QTIL    : nDOF x nM excitation vector in wavenumber domain
%   M       : integer wavenumber vector of the excitation ( exp(i*M*phi) )
%   F       : Frequency vector
%
% See also:

% Peter Fiala
% 2008 September

omvec = 2*pi*freqvec;
nDOF = size(K,1);
nOm = length(omvec);
nM = length(m);
ptil = zeros(nDOF, nM, nOm);
wb = waitbar(0, 'Computing response in the wavenumber domain');
for iOm = 1 : nOm
    om = omvec(iOm);
    for iM = 1 : nM
        ptil(:,iM,iOm) = (K + m(iM)^2*B - om^2*M) \ (-i*om*qtil(:,iM));
    end
    waitbar(iOm/nOm, wb);
end
close(wb);