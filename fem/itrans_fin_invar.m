function vz = itrans_fin_invar(vn, z)
%ITRANS_FIN_INVAR  Inverse modal expansion of response
%   Vz = ITRANS_FIN_INVAR(Vn, Z) Transforms the modal coordinates Vn of the
%   response into the spatial domain Vz.
%   Vn : modal coordinates of response (nDof x N+1 matrix)
%   Z  : M+1 element vector containing the values z0 z1 ... zM
%
% See also: trans_fin_invar, solve_finite_invariant

% Peter Fiala
% 2008

L = z(end)-z(1);
N = size(vn,2)-1;
A = sqrt(2.^(1-(0:N == 0)) / L);
nDof = size(vn,1);
vz = zeros(nDof,length(z));
for n = 0 : N
    vz = vz + vn(:,n+1) * A(n+1)*cos(n*pi/L*z);
end