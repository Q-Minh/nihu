function phat = solve_finite_invariant(K, M, A, csound, vz, z, N, freq)

%% Transform excitation to the wavenumber domain
vn = trans_fin_invar(vz, z, N);
qn = A * vn;

%% Solution in the wavenumber domain
ptil = zeros(size(vn));
L = z(end)-z(1);
om = 2*pi*freq;
for n = 0 : N
    ptil(:,n+1) = (K-(om^2-csound^2*(n*pi/L)^2)*M) \ (-i*om*qn(:,n+1));
end

%% Inverse transform into the spatial domain
phat = itrans_fin_invar(ptil, z);