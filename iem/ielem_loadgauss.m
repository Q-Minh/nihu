function G = ielem_loadgauss(type,ptype,P,param)
% Loads the gauss file for the given parameters

a = param(1);
b = param(2);

switch lower(ptype)
    case {1,'l','lag','lagrange'}
        gaussfile = sprintf('data/lag/gauss_ie%dm_%d',type-100,P);
    case {2,'j','jac','jacobi'}
        gaussfile = sprintf('data/jac%d%d/gauss_ie%dm_%d',a,b,type-100,P);
    otherwise
        error('ielem_mkc:invalidArg','Unknown polynomial family.');
end
gaussfile = fullfile(fileparts(mfilename('fullpath')), gaussfile);
% Load the required fields from the file
G = load(gaussfile,'w','N','dN','D','dD','dL','dMu');