function I = integpar(tree, k, C, acc, symm)
%INTEGPAR Compute integration parameters of a FM tree
%   I = integpar(tree, k, C, acc, symm) Computes the integration
%   parameters of a FM model. The integration parameters are computed for
%   each level of the cluster tree and involve:
% L      : truncation length
% S      : 3xQ Gaussian quadrature base points
% W      : Qx1 Gaussian quadrature weights
% Aup    : Interpolation matrices for upsampling
% Adn    : Interpolation matrices for downsampling
% Mz     : Qx1 z-symmetry matrix used for symmetric models
% Perm   : Qx16 Permutation matrices for symmetries
% Dindex : Translation distance indices for each interaction list pair
% Pindex : Permutation indices for each interaction list pair
% M      : Stored translation operators
%
% For symmetric structures, the parameters imDindex and imPindex are also
% stored.
% 
% Parameters:
%   tree : cluster tree obtained by clustertree
%   k    : Wave number
%   C    : accuracy parameter (between 2 and 4)
%
% See also: clustertree spherequad

% Peter Fiala
% 2009

%% Parameter check
if nargin < 5
    symm = 0;
end
if nargin < 4
    acc = 1;
end

%% Preallocation of integration parameter array
nL = size(tree,1); % number of levels in the tree (0 to d)
c = cell(nL,1);
I = struct('L', c,...   % expansion length
    'S', c,...          % quadrature points on the unit sphere
    'W', c,...          % quadrature weights
    'Aup', c,...        % Interpolation matrices
    'Adn', c, ...       % decimation matrices
    'Mz', c, ...        % Permutation matrix for z reflection
    'Dindex', c, ...    % Distance indices for each interlist element
    'imDindex', c, ...  % image distance indices for each interlist element
    'Perm', c, ...      % permutation matrix
    'Pindex', c, ...    % permutation indices for each interlist element
    'imPindex', c, ...  % image permutation indices for each interlist element
    'M', c);            % translation operators

%% Integration parameters for each level
mindepth = mininterdepth(tree, symm); % highest tranalation level
for l = nL : -1 : mindepth
    %% truncation length L
    x = sqrt(3)*abs(k)*tree(l).diameter;
    L = round(x+C*log10(x+pi));
    % TODO what is this for???
    if l < nL
        L = ceil(I(end).L + acc * (L - I(end).L));
    end
    I(l).L = 2*ceil(L/2); % ensure that L is even because of the quadrature
    
    %% quadrature over unit sphere
    [I(l).S, I(l).W, perm] = spherequad(I(l).L);
    I(l).Mz = perm(:,4);
    
    %% Interpolation matrices
    if l < nL
        [I(l+1).Aup I(l).Adn] = filter(I(l+1).L, I(l).L);
    end
    
    %% Unique distances with symmetry operators
    [D0, Dindex, I(l).Perm, Pindex] = Rsymmetry(tree(l).D, perm);
    Dindex = [0; Dindex]; %#ok<AGROW>
    Pindex = [0; Pindex]; %#ok<AGROW>
    I(l).Dindex = Dindex(tree(l).Dindex.'+1);
    I(l).Pindex = Pindex(tree(l).Dindex.'+1);
    if symm
        I(l).imDindex = Dindex(tree(l).imDindex.'+1);
        I(l).imPindex = Pindex(tree(l).imDindex.'+1);
    end
    
    %% Translation operators
    I(l).M = translation(k, D0, I(l).L, I(l).S);
end
end

%%
function [R0, Rindex, T0, Tindex] = Rsymmetry(R, perm)

R1 = abs(R);
ind1 = R1;
ind1(ind1 ~= 0) = ind1(ind1~=0) ./ R(ind1~=0);

T1 = repmat(perm(:,1), 1 ,size(R,1));
for j = 1 : 3
    a = find(ind1(:,j) == -1);
    T1(:,a) = T1(perm(:,j+1),a);
end
[R1, ~, n1] = unique(R1, 'rows');

%%
R2 = R1;
[R2(:,1:2), i] = sort(R2(:,1:2), 2);
ind2 = find(i(:,1) == 2);
T2 = repmat(perm(:,1), 1, size(R2,1));
for i = 1 : length(ind2)
    T2(:,ind2(i)) = T2(perm(:,5),ind2(i));
end
[R0, ~, n2] = unique(R2, 'rows');
Rindex = n2(n1);

%%
T0 = zeros(size(T2,1),length(Rindex));
for i = 1 : length(Rindex)
    T0(:,i) = T2(T1(:,i),n1(i));
end
[T0, ~, Tindex] = unique(T0.', 'rows');
T0 = T0.';
end

%%
function [Aup Adn] = filter(L1, L2)

%% Legendre functions used for interpolation
P1 = zeros(L1+1, L1+1, L1);
[z1 w1] = gaussquad(L1);
for m = 0 : L1
    P1(m+1,1:m+1,:) = legendre(m, z1, 'norm');
end
P2 = zeros(L2+1, L2+1, L2);
[z2 w2] = gaussquad(L2);
for m = 0 : L2
    P2(m+1,1:m+1,:) = legendre(m, z2, 'norm');
end

J = min(L1, L2);
Aup = zeros(J+1, L2, L1);
Adn = zeros(J+1, L1, L2);
for m = 0 : J
    iM = m+1;
    a = zeros(L2, L1);
    for l = abs(m) : J
        iL = l+1;
        a = a + squeeze(P2(iL,iM,:)) * squeeze(P1(iL,iM,:)).';
    end
    Aup(iM,:,:) = a * diag(w1);
    Adn(iM,:,:) = a.' * diag(w2);
end
end
