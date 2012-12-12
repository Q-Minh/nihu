function ptil = solve_infinite_invariant_direct(K, M, csound, qtil, kzvec, freqvec)

omvec = 2*pi*freqvec;
nDOF = size(K,1);
nOm = length(omvec);
nkz = length(kzvec);
ptil = zeros(nDOF, nkz, nOm);
wb = waitbar(0, 'Computing response in the wavenumber domain');
for iOm = 1 : nOm
    om = omvec(iOm);
    for ikz = 1 : nkz
        ptil(:,ikz,iOm) = (K - (om^2 - csound^2*kzvec(ikz)^2)*M) \ (-i*om*qtil(:,ikz));
    end
    waitbar(iOm/nOm, wb);
end
close(wb);