function gauss_gen_ie34m(P,ptype,alpha,beta)
%GAUSS_GEN_IE22M Generates element function data for two dimensional, Pth
%order infinite element based on jacobian polynomials with parameters
%(alpha,beta).

for c1 = 1:length(P)
    [ps, ws] = gaussquad(P(c1)+2);  % Gaussian quadrature in s direction
    [pt, wt] = gaussquad(2);    % Gaussian quadrature in t direction
    [pu, wu] = gaussquad(2);    % Gaussian quadrature in u direction
    [S, T, U] = ndgrid(ps,pt,pu); % Gauss points in the element
    xi = [S(:) T(:) U(:)].';    % Gauss points as a list
    W = ws*wt.';                % Gauss weights in the element
    W = W(:)*wu.';
    w = W(:);                   % Gauss weights as a list
     switch ptype
        case{'l','lag','lagrange'}
            gaussfile = sprintf('data/lag/gauss_ie34m_%d',P(c1));
            alpha = 0;
            beta = 0;
        case{'j','jac','jacobi'}
            gaussfile = sprintf('data/jac%d%d/gauss_ie34m_%d',alpha,beta,P(c1));
    end
    %Calculate element function values
    [N,dN,D,dD,dL,dMu] = ieshape34m(xi,P(c1),ptype,alpha,beta); %#ok<NOUSE>
    %Save matrices into the data file
    save(gaussfile,'w','N','dN','D','dD','dL','dMu');
end
