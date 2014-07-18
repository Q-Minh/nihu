clc
clear mex

make bemHG_const_bm.mex.cpp
make bemHG_const_bm_sp.mex.cpp
make recover_bm.mex.c

movefile(sprintf('*.%s', mexext), '../');

% make ibemB_const
