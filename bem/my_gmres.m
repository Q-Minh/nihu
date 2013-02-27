function [x, eps, times] = my_gmres(A, b, m, tol, Mright, varargin)
%MY_GMRES GMRES algorithm with right preconditioner
%   [x, eps, times] = my_gmres(A, b, m, tol, Mright, varargin) computes the
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
%   eps : Nit vector containing the total running time for each iteration
%
% MY_GMRES(..., 'quiet', true) suppresses textual output
%
% See also: sairight

% Peter Fiala
% 2013

%% extract options and check for -quiet
[opts, inner_args] = options(varargin{:});
if isfield(opts, 'quiet')
    quiet = opts.quiet;
else
    quiet = false;
end

%% allocating space for the algorithm
N = size(b,1);
h = zeros(m+1,m);   % Hessenberg matrix
v = zeros(N, m+1);  % Krylov subspace base vectors
eps = zeros(m, 1);  % residuals
times = zeros(m, 1);    % running time for each iteration

%% initialization phase
% TODO CHECK HAPPY CONDITION
x0 = zeros(N,1);
r0 = b - A(x0, inner_args{:});    % starting residual
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
    gmres_tic_id = tic();
    % Generate othonormal base of the Krylov subspace
    w = A(Mright * v(:,j), inner_args{:}); % new base vector
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
    times(j) = toc(gmres_tic_id);
    if ~quiet
        fprintf(1, '\tIteration: %d, residual: %e\n', j, eps(j));
    end
    % check convergence
    if eps(j) < tol
        break;
    end
end

eps = eps(1:j);
times = times(1:j);

%% build solution
y = R(1:end-1,:) \ g(1:end-1);
x = x0 + (Mright * v(:,1:j) * y);
end
