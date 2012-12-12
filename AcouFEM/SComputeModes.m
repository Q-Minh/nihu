function data = SComputeModes(data, nModes, nMaxIter)

ModalOpt.issym = 1;
ModalOpt.isreal = 1;
ModalOpt.maxit = nMaxIter;
ModalOpt.disp = 1;
ModalOpt.cholB = 0;

[data.Modes.Phi, data.Modes.Om] = fe_modes(data.Matrices.M, data.Matrices.K, nModes, ModalOpt);
data.Modes.DOF = data.Matrices.DOF;