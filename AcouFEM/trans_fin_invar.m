function vn = trans_fin_invar(vz, z, N)
%TRANS_FIN_INVAR Modal expansion of excitation
%   Vn = TRANS_FIN_INVAR(Vz, Z, N) Transforms the velocity excitation Vz
%   into a superposition of longitudinal modes.
%   Vz : Velocity excitation (nDof x M+1 matrix)
%   Z  : M+1 element vector containing the values z0 z1 ... zM
%   N  : The modes are computed within the index set n = 0 1 2 ... N

nDof = size(vz,1);
K = size(vz,2)-1;

L = z(end)-z(1);
dz = L / K;
q = (1:N)*pi/L;
qdz = q*dz;
vz0 = vz(:,1:end-1);
dv = diff(vz,1,2);

vn = zeros(nDof, N+1);
for k = 1 : K
    vn(:,1) = vn(:,1) + (vz0(:,k) + 0.5*dv(:,k)) * dz;
    vn(:,2:end) = vn(:,2:end) + vz0(:,k) * ((sin(qdz*k)-sin(qdz*(k-1)))./q) + ...
        dv(:,k) * ((cos(qdz*k)+qdz.*sin(qdz*k)-cos(qdz*(k-1)))/dz./q.^2);
end

A = sqrt(2.^(1-(0:N == 0)) / L);
vn = vn .* repmat(A, nDof, 1);