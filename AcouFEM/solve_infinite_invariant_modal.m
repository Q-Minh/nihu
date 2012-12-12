function [ptil, alphatil] = solve_infinite_invariant_modal(Phi, Om, csound, qtil, kzvec, freqvec)

omvec = 2*pi*freqvec;
nOm = length(omvec);
nMod = size(Om,1);
nDOF = size(Phi,1);
nkz = length(kzvec);
alphatil = zeros(nMod, nkz, nOm);
ptil = zeros(nDOF, nkz, nOm);
wb = waitbar(0, 'Computing response in the wavenumber domain');
Lambda = diag(Om.^2);
I = speye(size(Lambda));
PhiTqtil = Phi.' * qtil;
for iOm = 1 : nOm
    om = omvec(iOm);
    for ikz = 1 : nkz
        alphatil(:,ikz,iOm) = (Lambda - (om^2 - csound^2*kzvec(ikz)^2)*I) \ (-i*om*PhiTqtil(:,ikz));
    end
    ptil(:,:,iOm) = Phi * squeeze(alphatil(:,:,iOm));
    waitbar(iOm/nOm, wb);
end
close(wb);