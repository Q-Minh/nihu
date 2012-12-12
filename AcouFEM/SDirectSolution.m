function data = SDirectSolution(data)

K = data.Matrices.K;
M = data.Matrices.M;
switch data.ModelType
    case {'2DCart', '3DCart', '2DCyl'}
        qhat = zeros(size(data.Matrices.DOF));
        for iV = 1 : length(data.Excitation.NormalVelocity)
            qhat = qhat + data.Excitation.NormalVelocity(iV).A * data.Excitation.NormalVelocity(iV).vhat;
        end
        freq = data.Excitation.FrequencyVector;
        phat = solve_direct(K, M, qhat, freq);
        data.Response.phat = phat;
    case {'2.5DCartInf'}
        qtil = zeros(size(data.Matrices.DOF,1), length(data.Model.kz));
        for iV = 1 : length(data.Excitation.NormalVelocity)
            qtil = qtil + data.Excitation.NormalVelocity(iV).A * data.Excitation.NormalVelocity(iV).vtil;
        end
        freq = data.Excitation.FrequencyVector;
        csound = data.Model.Domain.Materials(1,4);
        ptil = solve_infinite_invariant_direct(K, M, csound, qtil, data.Model.kz, freq);
        data.Response.ptil = ptil;
    case {'2.5DCylFull'}
        B = data.Matrices.B;
        qtil = zeros(size(data.Matrices.DOF,1), length(data.Model.m));
        for iV = 1 : length(data.Excitation.NormalVelocity)
            qtil = qtil + data.Excitation.NormalVelocity(iV).A * data.Excitation.NormalVelocity(iV).vtil;
        end
        freq = data.Excitation.FrequencyVector;
        ptil = solve_full_cylindrical_direct(K, M, B, qtil, data.Model.m, freq);
        data.Response.ptil = ptil;
end