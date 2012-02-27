function [Ip, Iv] = fe_eval_interp(omega, Ap, Kp, Av, Kv)
% FE_EVAL_INTERP Evaluate interpolation matrices for sound pressure and 
%   particle velocity.
% IP = FE_EVAL_INTERP(OMEGA, AP, KP) returns the pressure interpolation
%   matrix IP for the given angular frequency OMEGA, magnitude matrix AP
%   and phase matrix KP. AP and KP can be determined for a given field
%   point mesh by using fe_invmap and fe_interp. 
% [IP, IV] = FE_EVAL_INTERP(OMEGA, AP, KP, AV, KV) returns the pressure and
%   particle velocity interpolation matrices IP and IV for the given
%   angular frequency OMEGA, pressure and velocity magnitude and phase
%   matrices AP, KP, AV, KV respectively. These matrices can be determined
%   for a given field point mesh by using fe_invmap and fe_interp.
% Final results for the sound pressure and particle velocity can be
% determined by the products Ip*p and Iv*p, where p is a column vector of
% the nodal pressure values.
%
% See also fe_invmap, fe_interp

% Peter Rucz, 2010.

% Check arguments
error(nargchk(3,5,nargin,'struct'));
error(nargoutchk(1,2,nargout,'struct'));

kAp = find(Ap);
[iAp, jAp] = find(Ap);
Ip = sparse(iAp,jAp,Ap(kAp).*exp(-1i*omega*Kp(kAp)),size(Ap,1),size(Ap,2));

if nargout > 1
    dim = size(Av,1) / size(Ap,1);
    Apv = reshape(repmat(Ap.',dim,1),size(Ap,2),dim*size(Ap,1)).';
    Kpv = reshape(repmat(Kp.',dim,1),size(Kp,2),dim*size(Kp,1)).';
    kAv = find(Av);
    [iAv, jAv] = find(Av);
    Iv = 1/(-1i*omega)*sparse(iAv,jAv, Av(kAv).*exp(-1i*omega*Kpv(kAv)) + ...
                         Apv(kAv).*(-1i*omega*Kv(kAv).*exp(-1i*omega*Kpv(kAv))),size(Av,1),size(Av,2)); 
end
end
