function [Phi, Om] = fe_modes(M, K, nModes, Opts)
%FE_MODES  Compute mode shapes and eigenfrequencies of a FE model
%   [PHI, OM] = FE_MODES(M, K, NMODES)
%   [PHI, OM] = FE_MODES(M, K, NMODES, OPTS)
%   Computes the lowest eigenfrequencies and corresponding mode shapes
%   of an acoustical FE model defined by its mass M and stiffness K
%   matrices.
%   PHI    : nDOF x nModes matrix containg the mode shapes
%   OM     : nModes element column vector containing the eigenfrequencies
%            in rad/s (2*pi*f)
%   NMODES : desired number of eigenfrequencies
%   OPTS   : eigenfunction computing options described in function EIGS
%
% See also: eigs, solve_modal, solve_infinite_invariant_modal,
% solve_full_cylindrical_modal

% Peter Fiala
% 2009

%% Default value for modal analysis options
if nargin < 4
    Opts.issym = 1;
    Opts.isreal = 1;
    Opts.maxit = 300;
    Opts.disp = 0;
    Opts.cholB = 0;
end

%% Eigenfunction computation
[Phi, Lam] = eigs(K, M, nModes, 'SM', Opts);

%% Postprocessing (sorting eigenfrequencies)
Lam = diag(Lam);
[Om, ind] = sort(real(sqrt(Lam)));
Phi = Phi(:,ind);
end
