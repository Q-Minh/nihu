function [x, eps] = my_gmres(A, b, m, tol, Mright, varargin)
%MY_GMRES GMRES algorithm with right preconditioner
%   [x, eps] = my_gmres(A, b, m, tol, Mright, varargin) computes the
%   iterative solution of the linear system  A x  = b with a right
%   preconditioner.
% Parameters:
%   A      : function handle that computes Ax
%   b      : right hand side vector
%   m      : max number of iterations
%   tol    : backward relative tolerance
%   Mright : Sparse Right preconditioner
% Output:
%   x   : Solution vector
%   eps : Nit vector containing the backward errors in each iteration step
%
% See also: sairight

% Peter Fiala
% 2009

%% allocating space for the algorithm
N = size(b,1);
h = zeros(m+1,m);   % Hessenberg matrix
v = zeros(N, m+1);  % Krylov subspace base vectors
eps = zeros(m, 1);  % residuals

%% initialization phase
% TODO CHECK HAPPY CONDITION
x0 = zeros(N,1);
r0 = b - A(x0, varargin{:});    % starting residual
beta = norm(r0);    % normalize
if beta < tol
    x = x0;
    eps = tol;
    return
end
v(:,1) = r0 / beta;

%% iteration phase
e1 = zeros(m+1,1);
e1(1) = 1;
for j = 1 : m
    % Generate othonormal base of the Krylov subspace
    w = A(Mright * v(:,j), varargin{:}); % new base vector
    % Gram-Schmidt orthogonalization
    ind = 1 : j;
    h(ind,j) = v(:,ind)' * w;
    for i = 1 : j
        w = w - v(:,i) * h(i,j);
    end
    % nomalization
    h(j+1,j) = norm(w);
    v(:,j+1) = w / h(j+1,j);
    % ls solution
    [Q, R] = qr(h(1:j+1,1:j));
    g = Q \ (beta * e1(1:j+1));
    eps(j) = abs(g(end))/beta;
    fprintf(1, '\tIteration: %d, residual: %e\n', j, eps(j));
    % check convergence
    if eps(j) < tol
        break;
    end
end

eps = eps(1:j);

%% build solution
y = R(1:end-1,:) \ g(1:end-1);
x = x0 + (Mright * v(:,1:j) * y);
end
