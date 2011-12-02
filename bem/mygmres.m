function [x,flag,relres,resvec] = mygmres(A,b,inner,tol,outer,M,varargin)
%MYGMRES GMRES iterative solver for FMBEM
%   [x,flag,relres,resvec] = MYGMRES(A,b,inner,tol,outer,M,varargin)
%   Slightly modified version of the Matlab strandard GMRES function
%   Main modification: Preconditioner M is used as premultiplier and not as
%   predivier:
%   Solve
%     MAx = Mb (current version)
%   instead of
%     M\(Ax) = M\b (oroginal GMRES version)
%
% See also: GMRES

% Peter Fiala
% 2009

n2b = norm(b);
n = size(b,1);
x = zeros(n,1);

% Set up for the method
flag = 1;
xmin = x;                        % Iterate which has minimal residual so far
jmin = 0;                        % "Inner" iteration at which xmin was computed
tolb = tol * n2b;                % Relative tolerance
if isa(A, 'function_handle')
    r = b - A(x,varargin{:});
else
    r = b - A*x;
end
normr = norm(r);                 % Norm of residual

if (normr <= tolb)               % Initial guess is a good enough solution
    flag = 0;
    relres = normr / n2b;
    resvec = normr;
    return
end

r = M * r;

normr = norm(r);
resvec = zeros(inner*outer+1,1);  % Preallocate vector for norm of residuals
resvec(1) = normr;                % resvec(1) = norm(b-A*x0)
normrmin = normr;                 % Norm of residual from xmin

%  Preallocate J to hold the Given's rotation constants.
J = zeros(2,inner);

for outiter = 1 : outer
    %  Construct u for Householder reflector.
    %  u = r + sign(r(1))*||r||*e1
    u = r;
    beta = scalarsign(r(1))*normr;
    u(1) = u(1) + beta;
    u = u / norm(u);

    U = zeros(n,inner);
    Psol = zeros(n,inner);
    R = zeros(n,inner);
    U(:,1) = u;

    %  Apply Householder projection to r.
    %  w = r - 2*u*u'*r;
    w = zeros(n,1); w(1) = -beta;

    for initer = 1 : inner
        %  Form P1*P2*P3...Pj*ej.
        %  v = Pj*ej = ej - 2*u*u'*ej
        v = -2*(u(initer)')*u;
        v(initer) = v(initer) + 1;
        for k = (initer-1):-1:1
            v = v - 2*U(:,k)*(U(:,k)'*v);
        end

        Psol(:,initer) = v;

        %  Apply preconditioner and A to v.
        if(isa(A, 'function_handle'))
            v = M * A(v,varargin{:});
        else
            v = M * A * v;
        end
        %  Form Pj*Pj-1*...P1*Av.
        for k = 1:initer
            v = v - 2*U(:,k)*(U(:,k)'*v);
        end

        %  Determine Pj+1.
        if (initer ~= length(v))
            %  Construct u for Householder reflector Pj+1.
            u = zeros(n,1);
            vhat = v(initer+1:end);
            alpha = norm(vhat);
            if (alpha ~= 0)
                alpha = scalarsign(vhat(1))*alpha;
                %  u = v(initer+1:end) +
                %        sign(v(initer+1))*||v(initer+1:end)||*e_{initer+1)
                u(initer+1:end) = vhat;
                u(initer+1) = u(initer+1) + alpha;
                u = u / norm(u);
                U(:,initer+1) = u;

                %  Apply Pj+1 to v.
                %  v = v - 2*u*(u'*v);
                v(initer+2:end) = 0;
                v(initer+1) = -alpha;
            end
        end

        %  Apply Given's rotations to the newly formed v.
        for colJ = 1:initer-1
            tmpv = v(colJ);
            v(colJ)   = conj(J(1,colJ))*v(colJ) + conj(J(2,colJ))*v(colJ+1);
            v(colJ+1) = -J(2,colJ)*tmpv + J(1,colJ)*v(colJ+1);
        end

        %  Compute Given's rotation Jm.
        if ~(initer==length(v))
            rho = norm(v(initer:initer+1));
            J(:,initer) = v(initer:initer+1)./rho;
            w(initer+1) = -J(2,initer).*w(initer);
            w(initer) = conj(J(1,initer)).*w(initer);
            v(initer) = rho;
            v(initer+1) = 0;
        end

        R(:,initer) = v;

        if initer < inner
            normr = abs(w(initer+1));
            resvec( (outiter-1)*inner+initer+1 ) = normr;
        end

        if normr <= normrmin
            normrmin = normr;
            jmin = initer;
        end

        fprintf(1, 'iteration %d - %d, residual: %x\n', outiter, initer, normr/n2b);

        if normr < tolb
            flag = 0;
            break
        end
    end         % ends innner loop
    
    y = R(1:jmin,:) \ w(1:jmin);
    additive = Psol*y;
    x = x + additive;
    xmin = x;
    if(isa(A, 'function_handle'))
        r = M * (b - A(x,varargin{:}));
    else
        r = M * (b - A*x);
    end

    normr = norm(r);

    resvec((outiter-1)*inner+initer+1) = normr;

    if normr <= normrmin
        xmin = x;
        normrmin = normr;
    end

    %  Test for stagnation on outer iterate.
    if flag~=2
        stagtest = zeros(n,1);
        ind = (x ~=0 );
        stagtest(ind) = additive(ind) ./ x(ind);
        if ( norm(stagtest,inf) < eps )
            flag = 3;
            break;
            %  No change in outer iterate.
        end
    end

    if normr < tolb
        flag = 0;
        break;
    end
end         % ends outer loop

% returned solution is that with minimum residual
if flag == 0
    relres = normr / n2b;
else
    x = xmin;
    relres = normrmin / n2b;
end

% truncate the zeros from resvec
if flag <= 1 || flag == 3
    resvec = resvec(1:(outiter-1)*inner+initer+1);
    indices = resvec==0;
    resvec = resvec(~indices);
else
    if initer == 0
        resvec = resvec(1:(outiter-1)*inner+1);
    else
        resvec = resvec(1:(outiter-1)*inner+initer);
    end
end

function sgn = scalarsign(d)
sgn = sign(d);
sgn(sgn == 0) = 1;