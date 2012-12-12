function data = SNewModel(data, type, varargin)

data = SInitData(data);
data.ModelType = type;
switch type
    case '2.5DCartFin'
        Zmax = varargin{1};
        Nz = varargin{2};
        data.Model.zVector = linspace(0, Zmax, Nz+1);
    case '2.5DCartInf'
        Zmin = varargin{1};
        Zmax = varargin{2};
        Nz = varargin{3};
        data.Model.zVector = linspace(Zmin, Zmax, Nz+1);
    case '2.5DCylFin'
        Phimax = varargin{1};
        nPhi = varargin{2};
        data.Model.phiVector = linspace(0, Phimax, nPhi+1);
    case '2.5DCylFull'
        nPhi = varargin{1};
        data.Model.phiVector = linspace(0, 2*pi, nPhi+1);
end