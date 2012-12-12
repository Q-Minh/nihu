function data = STransformSolution(data)

switch data.ModelType
    case {'2.5DCartInf'}
        [phat] = itrans_inf_invar(data.Response.ptil, data.Model.kz);
        s = size(phat);
        data.Response.phat = reshape(phat, prod(s(1:2)), s(3));
    case {'2.5DCylFull'}
        [phat] = itrans_full_cyl(data.Response.ptil, data.Model.m);
        s = size(phat);
        data.Response.phat = reshape(phat, prod(s(1:2)), s(3));
    otherwise
end