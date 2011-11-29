function f = ispht(Flm, Plmj)
%ISPHT   Inverse spherical harmonic transform
%   F = ISPHT(Flm, Plmx) Computes the inverse spherical harmonic transform
%   of the function Flm.
%   Flm  : 3D matrix containing the amplitudes of spherical harmonics
%          Y_l^m. The first dimension is l = 0:L, the second is m = -L:L,
%          and the third dimension is q = 1:Q. The transformation is
%          performed independently for each value of q.
%   Plmj : 3D matrix of associated Legendre functions P_l^m(x_j). The
%          first dimension is l = 0:L, the second is m = 0:l and the third
%          is j = 1:L, and x_j denotes the Gaussian quadrature points
%          between -1 and +1.
%   F    : 3D matrix containing the inverse transformed function samples.
%          The first dimension is j = 1:L, the second dimension is p = 1:2L
%          and the third is q=1:Q. The result is returned in the sample
%          points F(theta_j,phi_p) where theta_t is acos(x_j) and phi_p =
%          (p-1)/L*pi.
%
% See also: spht, sphresamp

L2 = size(Flm,1)-1; % length of spectrum
Q = size(Flm, 3);
L1 = size(Plmj,1)-1;

%% inverse transfrom along the 1. (theta) direction
Phi_m = zeros(L1, 2*L1, Q);
l = (0:L2).';
m = -L2 : L2-1;
m1 = L1+1+m;
m2 = L2+1+m;
plmj = Plmj(l+1,abs(m)+1,:);
ilflm = repmat(i.^l,[1,length(m),Q]) .* Flm(l+1,m2,:);
for j = 1 : L1
    Phi_m(j,m1,:) = squeeze(sum(repmat(plmj(:,:,j), [1 1 Q]) .* ilflm, 1));
end

%% inverse transfrom along the 2. (phi) direction
f = 2*L1*ifft(fftshift(Phi_m,2),[],2);