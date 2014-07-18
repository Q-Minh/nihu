function varargout = model_system(varargin)
%   [M K DOF] = [indsys1, indsys2, ...]
%   [M K C DOF] = [indsys1, indsys2,...]

I = [];
J = [];
M = [];
K = [];
C = [];
DOF = [];
l = 0;

for iSys = 1:nargin
    I = [I; varargin{iSys}{1}+l];   %#ok<AGROW>
    J = [J; varargin{iSys}{2}+l];   %#ok<AGROW>    
    M = [M; varargin{iSys}{3}];     %#ok<AGROW>
    K = [K; varargin{iSys}{4}];     %#ok<AGROW>
    if length(varargin{iSys}) == 5;
        DOF = [DOF; varargin{iSys}{5}];     %#ok<AGROW>
    else
        C = [C; varargin{iSys}{5}];         %#ok<AGROW>
        DOF = [DOF; varargin{iSys}{6}];     %#ok<AGROW>
    end
    l = length(DOF);
end

% Unique the DOF (also does the sorting)
[DOF, uind, rind] = unique(DOF);            %#ok<ASGLU>
I = rind(I);        % Re-index row indices
J = rind(J);        % Re-index column indices
l = length(DOF);

% Assemble sparse matrices
M = sparse(I, J, M, l, l);
K = sparse(I, J, K, l, l);

varargout{1} = M;
varargout{2} = K;
if nargout == 3
    varargout{3} = DOF;
elseif nargout == 4
    C = sparse(I, J, C, l, l);
    varargout{3} = C;
    varargout{4} = DOF;
end

end % of function model_system

