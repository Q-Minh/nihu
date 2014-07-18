function f2 = sphresamp(f, L1, L2, A)
%SPHRESAMP   Resampling a function over the unit sphere
%   F2 = SPHRESAMP(F, L1, L2, A)  resamples the function F defined over the
%   unit sphere. The original function F is defined in the integration base
%   points contained in the structure I1, and the resampled function F2
%   will be written over the new integration base points given in I2. The
%   resampling algorithm is based on a spherical harmonic transform and
%   zero padding for upsampling and truncation for downsampling.
% Parameters:
%   F      : m x n complex matrix, each column contains the samples of a
%            function over the unit sphere. The base points are contained
%            in I1.S. The resampling is performed simultaneously for each
%            column of F.
%   L1, L2 : Integration parameter structures obtained by INTEGPAR
%   F2     : M x n complex matrix, each column contains the samples of the
%            resampled function. The corresponding base points are I2.S.

%% shortcut
if L1 == L2
    f2 = f;
    return;
end

%% Parameters
Q = size(f,2);

%% FFT along the 2. (phi) direction
Phi_m = 1/(2*L1) * fftshift(fft(reshape(f, L1, 2*L1, Q), [], 2), 2);

%% Filtering along the 1. (theta / z) dimension
J = min(L1, L2);
Phi = zeros(L2, 2*L2, Q);
for m = -J : J-1
    Phi(:,m+L2+1,:) = ...
        squeeze(squeeze(A(abs(m)+1,:,:))) * squeeze(Phi_m(:,m+L1+1,:));
end

%% Inverse FFT along the 2. (phi) direction
f2 = reshape(2*L2*ifft(fftshift(Phi,2),[],2), 2*L2^2, Q);
end