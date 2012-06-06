function Hp = gmres_iter_neumann(p, Hnf, rr, rs, ns, w, gind, tr, intdata, k, fs, fr, symm)
%GMRES_ITER_NEUMANN One GMRES iteration of a Neumann problem
%   Hp = gmres_iter_neumann(p, Hnf, rr, rs, ns, w, gind, tr, intdata, k,
%   fs, fr, symm) computes one iteration of a GMRES algorithm used to
%   compute the discretized integral operator Hp.
%   The sources are the Gaussian integration points of the mesh.
%   The receivers are the DOF (element center) of the bem mesh.
% Parameters:
%   p       : Nx1 last iteration result
%   Hnf     : NxN near field sparse matrix
%   rr      : 3xN matrix, xyz coordinates of receivers.
%   rs      : 3xM matrix, xyz coordinates of sources.
%   ns      : 3xM matrix, xyz coordinates of source normals
%   w       : Mx1 vector of source weights (Gaussian quadrature)
%   gind    : Mx1 Index vector containing the element index of Gaussian points.
%   rtree   : relative cluster tree obtained by reltree
%   intdata : integration parameter structure built by integpar
%   k       : scalar real wave number
%   fs      : father cluster indices of source elements returned by
%             clustertree
%   fr      : father cluster indices of receiver elements returned by
%             clustertree
%   symm    : symmetry parameter (-1, 0 or +1) (antisymm, no symm, symm)
%
% See also: clustertree reltree integpar mpcont_Hp

% Peter Fiala
% 2009

%% parameter check
if nargin < 13
    symm = 0;
end

%% Surface integral
Hp = -.5*p + Hnf*p + ...
    mpcont_Hp(rr, rs, ns, p(gind).*w, tr, intdata, k, fs(gind), fr, symm);
end
