function data = SComputeMK(data)

switch data.ModelType
    case {'3DCart', '2DCart', '2.5DCartInf', '2.5DCartFin'}
        [data.Matrices.M, data.Matrices.K, data.Matrices.DOF] = model_mk(data.Model.Domain);
    case {'2DCyl', '2.5DCylFull', '2.5DCylFin'}
        [data.Matrices.M, data.Matrices.K, data.Matrices.B, data.Matrices.DOF] = model_mk_cyl(data.Model.Domain);
end