function Flm = spht(f, w, Plmx, L2)
%SPHT   Spherical harmonic transform
%   Flm = SPHT(F, W, Plmj, L2) Computes the spherical harmonic transform
%   of the function F.
%   F    : 3D matrix containing the input function samples. The first
%          dimension is j = 1:L, the second dimension is p = 1:2L
%          and the third is q=1:Q. The function is sampled in the points
%          F(theta_j,phi_p) where theta_t is acos(x_j) and x_j denotes the
%          Gaussian quadrature points between -1 and +1, and phi_p =
%          (p-1)/L*pi. The transformation is performed independently for
%          each value of q.
%   W    : L-point Gaussian weights between -1 and +1.
%   Plmj : 3D matrix of associated Legendre functions P_l^m(x_j). The
%          first dimension is l = 0:L, the second is m = 0:l and the third
%          is j = 1:L.
%   L2   : Number of harmonics to determine, the default value is L2 + L.
%          If L2 > L then zero padding is used, if L2 < L, truncation of the
%          spectrum is applied.
%   Flm  : 3D matrix containing the amplitudes of spherical harmonics
%          Y_l^m. The first dimension is l = 0:L2, the second is m = -L2:L2,
%          the third dimension is q=1:Q.
%
% See also: ispht, sphresamp

L = size(f,1);
Q = size(f, 3);

if nargin < 4
    L2 = L;
end

%% FFT along the 2. (phi) direction
Phi_m = 1/(2*L) * fftshift(fft(f, [], 2), 2);
Phi_m(:,end+1,:) = Phi_m(:,1,:);
clear f;

%% transform along the 1. (theta) direction
% l = 0 case is handled separately
Flm = zeros(L2+1, 2*L2+1, Q);
Flm(1,L2+1,:) = sqrt(2)/2 * (w.' * squeeze(Phi_m(:,L+1,:)));
% l > 0 case
S = min([L, L2]);
for l = 1 : S
    m = -l : +l;
    Flm(l+1,m+L2+1,:) = i^(-l) * ...
        squeeze(sum(...
        repmat(repmat(w,1,length(m)) .* squeeze(Plmx(l+1,abs(m)+1,:)).', [1 1 Q]) .* ...
        Phi_m(:,m+L+1,:),1));
end