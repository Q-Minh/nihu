function gauss_gen_ie33mj(P,ptype,alpha,beta)
%GAUSS_GEN_IE22M Generates element function data for two dimensional, Pth
%order infinite element based on jacobian polynomials with parameters
%(alpha,beta).

for c1 = 1:length(P)
    [ps, ws] = gaussquad(P(c1)+2);  % Gaussian quadrature in s direction
    [ptu, wtu] = dunavant_rule(2);    % Gaussian quadrature in t direction
    %Correct dunavant rules to match parent domain
    ptu(1,:) = -2*ptu(1,:)+1;
    ptu(2,:) =  2*ptu(2,:)-1;
    ntu = size(ptu,2);
    wtu      =  2*wtu;
    ps = reshape(repmat(ps.',ntu,1),(P(c1)+2)*ntu,1);
    ws = reshape(repmat(ws.',ntu,1),(P(c1)+2)*ntu,1);
    ptu = repmat(ptu.',P(c1)+2,1);
    wtu = repmat(wtu.',P(c1)+2,1);
    xi = [ps(:) ptu(:,1) ptu(:,2)].';    % Gauss points as a list
    w = ws.*wtu;                % Gauss weights in the element as a list
    switch ptype
        case{'l','lag','lagrange'}
            gaussfile = sprintf('data/lag/gauss_ie33m_%d',P(c1));
            alpha = 0;
            beta = 0;
        case{'j','jac','jacobi'}
            gaussfile = sprintf('data/jac%d%d/gauss_ie33m_%d',alpha,beta,P(c1));
    end
    %Calculate element function values
    [N,dN,D,dD,dL,dMu] = ieshape33m(xi,P(c1),ptype,alpha,beta); %#ok<NOUSE>#
    %Save matrices into the data file
    save(sprintf(gaussfile,'w','N','dN','D','dD','dL','dMu'));
end
