function ptil = solve_full_cylindrical_modal(Phi, Om, Gamma, qtil, m, freqvec)

omvec = 2*pi*freqvec;
nOm = length(omvec);
nMod = size(Om,1);
nDOF = size(Phi,1);
nM = length(m);
alphatil = zeros(nMod, nM, nOm);
ptil = zeros(nDOF, nM, nOm);
wb = waitbar(0, 'Computing response in the wavenumber domain');
Lambda = diag(Om.^2);
I = speye(size(Lambda));
PhiTqtil = Phi.' * qtil;
for iOm = 1 : nOm
    om = omvec(iOm);
    for iM = 1 : nM
        alphatil(:,iM,iOm) = (Lambda + m(iM)^2 * Gamma - om^2*I) \ (-i*om*PhiTqtil(:,iM));
    end
    ptil(:,:,iOm) = Phi * squeeze(alphatil(:,:,iOm));
    waitbar(iOm/nOm, wb);
end
close(wb);