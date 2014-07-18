function [U, S, delta] = excdecomp(b, eps, U, S, V)
%EXCDECOMP Decomposition of right hand side matrix for iterative solver
%   [U, S, delta] = excdecomp(b, eps)
%   [U, S, delta] = excdecomp(b, eps, U, S, V)
%   excdecomp decomposes the right-hand side matrix b of a linear
%   system A x = b into the form b = U * S and returns the tolerances of
%   the reduced system A y = U. The decomposition is based on a singular
%   value decomposition of b. If the svd has already been performed, the
%   call excdecomp(b, eps, U, S, V) is more efficient.
% Input parameters:
%   b       : n x m right-hand side of the linear system A x = b
%   eps     : 1 x m vector containing the prescribed tolerances for each
%             right-hand side vector of the system
%   U, S, V : svd decomposition matrices of the right-hand side vector b
% Output parameters:
%   U, S    : New decomposition of b = U * S, where U is an n x q matrix
%             and S is q x m.
%   delta   : 1 x q vector containing the tolerances of the reduced system
%             A * y = U

% Peter Fiala
% 2009

%% Singular value decomposition of B if needed
if nargin < 3
    [U, S, V] = svd(b, 0);
end
sigma = diag(S);
S = S * V';
clear V;

%% determine reduced base
beta = .4; % accuracy parameter, Heuristic choice (lang03a)
p = size(b,2);
normb = zeros(1,p);
for i = 1 : p
    normb(i) = norm(b(:,i));
end
normb = normb * eps;
q = max(find(sigma <= beta * min(normb), 1, 'first') - 1, 1);
U = U(:,1:q);
S = S(1:q,:);

%% determine tolerance for the reduced base
alpha = .4; % accuracy parameter, Heuristic choice (lang03a)
delta = zeros(1,q);
for j = 1 : q
    delta(j) = alpha * min(normb ./ abs(S(j,:)));
end
% drop tolerances larger than one
ind = find(delta < 1);
delta = delta(ind);
U = U(:,ind);
S = S(ind,:);
end
