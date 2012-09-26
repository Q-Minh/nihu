clc
clear mex

make bemHG_const_bm.mex.c
make bemHG_const_bm_sp.mex.c
make recover_bm.mex.c

movefile(sprintf('*.%s', mexext), '../');

% make ibemB_const