function M = translation(k, Dvec, L, s)
%TRANSLATION  Compute translation between multipoles
%   M = TRANSLATION(K, D, L, S)
%   K : wave number (scalar)
%   D : nD x 3 matrix (distance vectors in each row)
%   L : truncation parameter
%   S : 3 x nG matrix (Quadrature base points on unit sphere)

%% WARNING
% The original formula has been designed for exp(ikr)/r,
% the current formula is designed for exp(-ikr)/r
k = -k;

%% compute spherical Hankel functions
% find unique scalar distances
[D, ~, Dind] = unique(round(1e10*sqrt(dot(Dvec, Dvec, 2)))/1e10);
[Dir, ~, Dirind] = unique(round(1e10*(Dvec ./ repmat(D(Dind),1,3)))/1e10, 'rows');
h = sphankel((0:L), k*D);
%% compute translation matrix
x = Dir * s;
M = zeros(length(Dirind),size(x,2));
% initialize Legendre polynomials
P0 = ones(size(x));
P1 = x;
for l = 0:L
    % compute Legendre polynomials
    switch l
        case 0
            P = P0;
        case 1
            P = P1;
        otherwise
            % recursive evaluation of P_l
            P = ((2*l-1)*(x.*P1)-(l-1)*P0)/l;
            P0 = P1;
            P1 = P;
    end
    q = ((2*l+1)*(1i^l))*h(Dind,l+1);
    M = M + repmat(q,1,size(P,2)) .* P(Dirind,:);
    progbar(0, L, l);
end
M = M*(1i*k)/(4*pi)/(4*pi);
if isreal(k) && k < 0
    M = -M;
end

M = M.';

end

%% Spherical hankel functions
function h = sphankel(nu, z)

h = repmat(sqrt(pi./(2*z)),1,length(nu)) .* ...
    besselh(repmat(nu+1/2, length(z),1), repmat(z, 1, length(nu)));
end