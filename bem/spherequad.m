function [S, W, perm, model] = spherequad(p)
%SPHEREQUAD   Gaussian quadrature over the unit sphere
%   [S, W, x, w, model] = SPHEREQUAD(P)   returns Gaussian quadrature data
%   over the unit sphere. The quadrature scheme
%   Parameters:
%   P     : quadrature density parameter
%   S     : quadrature pase points (3xN) matrix, where N=(P*2P) is the
%           number of quadrature points, and each column gives the xyz
%           values of a quadrature points
%   W     : N x 1 column vector of quadrature weights
%   model : mesh of the unit circle for display purposes
%
% See also: gaussquad

% Peter Fiala
% 2009

%% Parameter check
pp = 2*ceil(p/2);
if pp ~= p
    warning('MATLAB:paramAmbiguous', 'The input parameter p has been changed to the nearest larger even integer.');
    p = pp;
end

%% Linear quadratures
phi = (0:2*p-1).'*pi/p; % uniform in phi direction
[z, w] = gaussquad(p);  % Gaussian in theta direction
theta = acos(z);

%% 2D quadrature
% quadrature base points
[Phi, Theta] = meshgrid(phi, theta);
Theta = Theta(:);
Phi = Phi(:);
S = [cos(Phi).*sin(Theta), sin(Phi).*sin(Theta), cos(Theta)].';
% quadrature weights
W = repmat(w, 2*p, 1) * pi/p;

%% Permutations
Iphi = 1 : 2*p;
Mxphi = mod(p : -1 : -p+1, 2*p)+1;
Myphi = mod(0 : -1 : -2*p+1, 2*p)+1;
Mxyphi = mod(p/2 - (0 : 2*p-1), 2*p)+1;
Itheta = (1 : p).';
Mztheta = (p : -1 : 1).';

perm = zeros(2*p^2,5);
perm(:,1) = perm2perm(Itheta, Iphi);    % identity
perm(:,2) = perm2perm(Itheta, Mxphi);   % Mx
perm(:,3) = perm2perm(Itheta, Myphi);   % My
perm(:,4) = perm2perm(Mztheta, Iphi);   % Mz
perm(:,5) = perm2perm(Itheta, Mxyphi);  % Mxy

%% Construct sphere geometry for plotting
if nargout == 4
    model = create_slab([1,1], [p-1,2*p-1]);
    model.Nodes(:,2:4) = S.';
    model.Elements(end+1:end+p-1,1:8) = [
        (p-1)*(2*p-1)+(1:p-1).', repmat([24 1 1],p-1,1),...
        repmat((1:p-1).',1,4) + repmat([(2*p-1)*p 0 1 (2*p-1)*p+1],p-1,1)
        ];
end
end

%%
function perm = perm2perm(theta, phi)

p = length(theta);
perm = reshape(1:2*p^2, p, 2*p);
perm = perm(theta, phi);
perm = perm(:);
end
