function [Ap, Kp, Av, Kv] = fe_interp(mesh, El, Xi)
%FE_INTERP Interpolation matrix based on inverse mapping
%  [Ap, Kp] = FE_INTERP(MESH, EL, XI) or
%  [Ap, Kp, Av, Kv] = FE_INTERP(MESH, EL, XI) returns NxM sparse 
%  interpolation matrices Ap and Av that relate the pressure and velocity 
%  in the N field points POINTS to the M DOF of a NiHu mesh. Matrices Kp
%  and Kv relate the phase terms for the pressure and particle velocity in
%  similar way. The complete interpolation matrices can be evaluated by the
%  function fe_eval_interp.
%
%See also fe_eval_interp.

% Peter Fiala, Peter Rucz
% Last modified: 16.04.2010.

%% Argument check and default parameters
error(nargchk(3,4,nargin,'struct'));
error(nargoutchk(1,4,nargout,'struct'));

%% Initialization
nNode = size(mesh.Nodes,1);         % number of nodes in the mesh
nPoint = size(El,1);                % number of field points
points = (1:nPoint).';
dim = floor(mod(mesh.Elements(1,2),100)/10); % dimensionality of the problem
% count all nodes
elements = drop_IDs(mesh);
ElOk = El(El > 0);                  % the non-zero element indices
% Check for infinite elements
ElInf = ElOk(elements(ElOk,2) > 100);
if isempty(ElInf)
    nodenum = sum(mod(elements(ElOk,2),10));
else
    ieP = mesh.Properties(elements(ElInf,4),3);
    nodenum = sum(mod(elements(ElOk,2),10)) + sum(mod(elements(ElInf,2),10).*(ieP-1));
end

%% Preallocate indices for sparse output matrices
rows_p = zeros(nodenum,1);
cols_p = zeros(nodenum,1);
Ap_v = zeros(nodenum,1);
Kp_v = zeros(nodenum,1);
% Allocate velocity matrices
if nargout > 2
    rows_v = zeros(nodenum*dim,1);
    cols_v = zeros(nodenum*dim,1);
    Av_v = zeros(nodenum*dim,1);
    Kv_v = zeros(nodenum*dim,1);
end
% Points are not ready
ready = false(nPoint,1);

%% If only the pressure information is needed
if nargout < 3
    ind_p = 0;
    for iPoint = 1:nPoint
        progbar(1,nPoint,iPoint);
        % If the point is already taken care for
        if(ready(iPoint) || El(iPoint) == 0)
            continue
        end
        iElem = (El == El(iPoint));
        type = elements(El(iPoint),2);
        xi = Xi(iElem,:);
        c = mesh.Materials(elements(El(iPoint), 3),4);
        switch type
            % Finite elements
            case {12, 23, 24, 34, 36, 38}
                nodes = elements(El(iPoint),5:(4+mod(type,10)));
                Ap_r = shapefun(xi,type);
                Kp_r = zeros(size(Ap_r));
            % Infinite elements
            case {111, 122, 133, 134}
                % Get properties of the element
                P = mesh.Properties(elements(El(iPoint),4),3);
                % Get element nodes
                nodes = elements(El(iPoint),[5:(mod(type,10)+4),...
                    5+mod(type,10)*2:4+(mod(type,10)*(P+1))]);
                % Get mapping coordinates
                map = mesh.Nodes(elements(El(iPoint),5:(mod(type,10)*2+4)),2:(dim+1));
                % Call postprocessing
                [Ap_r, Kp_r] = ie_postproc(type,map,P,xi);
        end
        % Create matrix indices
        ind_l = sum(iElem) * numel(nodes);
        Ap_v(ind_p+1:ind_p+ind_l) = Ap_r(:);
        Kp_v(ind_p+1:ind_p+ind_l) = Kp_r(:) / c;
        rows_p(ind_p+1:ind_p+ind_l) = repmat(points(iElem),numel(nodes),1);
        cols_p(ind_p+1:ind_p+ind_l) = reshape(repmat(nodes(:).',sum(iElem),1),ind_l,1);
        % Increase current index number
        ind_p = ind_p + ind_l;
        ready(iElem) = true;
    end
    Ap = sparse(rows_p,cols_p,Ap_v,nPoint,nNode);
    Kp = sparse(rows_p,cols_p,Kp_v,nPoint,nNode);
% If velocity information is also required
else
    ind_p = 0;
    ind_v = 0;
    % Check if point is already ready
    for iPoint = 1:nPoint
        progbar(1,nPoint,iPoint);
        if(ready(iPoint) || El(iPoint) == 0)
            continue
        end
        iElem = (El == El(iPoint));
        type = elements(El(iPoint),2);
        xi = Xi(iElem,:);
        c = mesh.Materials(elements(El(iPoint), 3),4);
        rho = mesh.Materials(elements(El(iPoint), 3),3);
        switch type
            case {12, 23, 24, 34, 36, 38}
                nodes = elements(El(iPoint),5:(4+mod(type,10)));
                xnod = mesh.Nodes(nodes,2:(dim+1));
                [N, dNxi] = shapefun(xi, type);
                Ap_r = N;
                Kp_r = zeros(size(Ap_r));
                Av_r = zeros(dim*size(Ap_r,1),numel(nodes));
                Kv_r = zeros(size(Av_r));
                for c1 = 1:size(xi,1)
                    Jac = squeeze(dNxi(c1,:,:)).';
                    Av_r(dim*(c1-1)+(1:dim), :) = (Jac * xnod) \ Jac;
                end
            case {111 122 133 134}
                 % Get properties of the element
                P = mesh.Properties(elements(El(iPoint),4),3);
                % Get element nodes
                nodes = elements(El(iPoint),[5:(mod(type,10)+4),...
                    5+mod(type,10)*2:4+(mod(type,10)*(P+1))]);
                % Get mapping coordinates
                map = mesh.Nodes(elements(El(iPoint),5:(mod(type,10)*2+4)),2:(dim+1));
                % Call postprocessing
                [Ap_r, Kp_r, Av_r, Kv_r] = ie_postproc(type,map,P,xi);
        end
        % Create matrix indices
        ind_l = sum(iElem) * numel(nodes);
        Ap_v(ind_p+1:ind_p+ind_l) = Ap_r(:);
        Kp_v(ind_p+1:ind_p+ind_l) = Kp_r(:) / c;
        Av_v(ind_v+1:ind_v+dim*ind_l) = Av_r(:) / rho;
        Kv_v(ind_v+1:ind_v+dim*ind_l) = Kv_r(:) / c / rho;
        rows_p(ind_p+1:ind_p+ind_l) = repmat(points(iElem),numel(nodes),1);
        rows_v(ind_v+1:ind_v+ind_l*dim) = repmat(reshape((repmat((points(iElem)-1)*dim,1,dim) + repmat((1:dim),sum(iElem),1)).',dim*sum(iElem),1),numel(nodes),1);
        cols_p(ind_p+1:ind_p+ind_l) = reshape(repmat(nodes(:).',sum(iElem),1),ind_l,1);
        cols_v(ind_v+1:ind_v+ind_l*dim) = reshape(repmat(nodes(:).',dim*sum(iElem),1),ind_l*dim,1);
        % Increase current index numbers
        ind_p = ind_p + ind_l;
        ind_v = ind_v + dim*ind_l;
        ready(iElem) = true;
    end
    Ap = sparse(rows_p,cols_p,Ap_v,nPoint,nNode);
    Kp = sparse(rows_p,cols_p,Kp_v,nPoint,nNode);
    Av = sparse(rows_v,cols_v,Av_v,nPoint*dim,nNode);
    Kv = sparse(rows_v,cols_v,Kv_v,nPoint*dim,nNode);
end
end