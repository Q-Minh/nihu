function data = SInitData(varargin)
%INITDATA
%   DATA = INITDATA
%   DATA = INITDATA(DATA)

if nargin == 1
    data = varargin{1};
end

data.ModelType = [];
data.Model = [];
data.Matrices = [];
data.Modes = [];
data.Excitation = [];
data.Response = [];