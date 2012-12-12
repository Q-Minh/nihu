function data = SComputeExcitationMatrix(data)

boundary = data.Model.Boundary;
NV = data.Excitation.NormalVelocity;
nV = length(NV);
for iV = 1 : nV
    bou = boundary;
    [nodeind, elemind] = mesh_select(bou, NV(iV).ElemSelector, 'ind');
    bou.Elements = bou.Elements(elemind,:);
    switch data.ModelType
        case {'3DCart', '2DCart', '2.5DCartFin', '2.5DCartInf'}
            NV(iV).A = model_a(bou);
        case {'2DCyl', '2.5DCylFull', '2.5DCylFin'}
            NV(iV).A = model_a_cyl(bou);
    end
end
data.Excitation.NormalVelocity  = NV;