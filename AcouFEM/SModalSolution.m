function data = SModalSolution(data)

Modes = data.Modes;
Excitation = data.Excitation;
switch data.ModelType
    case {'2DCart', '3DCart', '2DCyl'}
        qhat = zeros(size(data.Matrices.DOF));
        for iV = 1 : length(Excitation.NormalVelocity)
            qhat = qhat + Excitation.NormalVelocity(iV).A * Excitation.NormalVelocity(iV).vhat;
        end
        phat = solve_modal(Modes.Phi, Modes.Om, qhat, Excitation.FrequencyVector);
        data.Response.phat = phat;
    case {'2.5DCartInf'}
        qtil = zeros(size(data.Matrices.DOF,1), length(data.Model.kz));
        for iV = 1 : length(data.Excitation.NormalVelocity)
            qtil = qtil + data.Excitation.NormalVelocity(iV).A * data.Excitation.NormalVelocity(iV).vtil;
        end
        freq = data.Excitation.FrequencyVector;
        csound = data.Model.Domain.Materials(1,4);
        ptil = solve_infinite_invariant_modal(Modes.Phi, Modes.Om, csound, qtil, data.Model.kz, freq);
        data.Response.ptil = ptil;
    case {'2.5DCylFull'}
        B = data.Matrices.B;
        Gamma = data.Modes.Phi.' * B * data.Modes.Phi;
        qtil = zeros(size(data.Matrices.DOF,1), length(data.Model.m));
        for iV = 1 : length(data.Excitation.NormalVelocity)
            qtil = qtil + data.Excitation.NormalVelocity(iV).A * data.Excitation.NormalVelocity(iV).vtil;
        end
        freq = data.Excitation.FrequencyVector;
        ptil = solve_full_cylindrical_modal(Modes.Phi, Modes.Om, Gamma, qtil, data.Model.m, freq);
        data.Response.ptil = ptil;
end