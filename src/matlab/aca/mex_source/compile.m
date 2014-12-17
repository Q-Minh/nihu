eigendir = '/D/Peti/research/toolbox/eigen';
opts = strcat('-I', eigendir);

mex('M2L.cpp', opts, '-output', 'compute_M2L');
mex('M2M.cpp', opts, '-output', 'compute_M2M');
mex('L2L.cpp', opts, '-output', 'compute_L2L');

movefile *.mexa* ../