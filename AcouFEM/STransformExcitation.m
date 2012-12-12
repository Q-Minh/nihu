function data = STransformExcitation(data)

switch data.ModelType
    case {'2.5DCartInf'}
        zVector = data.Model.zVector;
        NV = data.Excitation.NormalVelocity;
        nV = length(NV);
        for iV = 1 : nV
            [NV(iV).vtil kz] = trans_inf_invar(reshape(NV(iV).vhat, size(NV(iV).vhat,1)/length(zVector), length(zVector)), zVector);
        end
        data.Excitation.NormalVelocity  = NV;
        data.Model.kz = kz;
    case {'2.5DCylFull'}
        phiVector = data.Model.phiVector;
        NV = data.Excitation.NormalVelocity;
        nV = length(NV);
        for iV = 1 : nV
            [NV(iV).vtil m] = trans_full_cyl(reshape(NV(iV).vhat, size(NV(iV).vhat,1)/(length(phiVector)-1), length(phiVector)-1), length(phiVector)-1);
        end
        data.Excitation.NormalVelocity  = NV;
        data.Model.m = m;
    otherwise
end