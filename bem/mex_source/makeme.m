clear mex

% conv bem
make bemHG_const.mex.cpp
make bemHG_const_sp.mex.cpp
make bemHG_lin.mex.c

% make ibemB_lin.mex.c
% make ibemD_const.mex.c
% make ibemD_lin.mex.c

% conv bem with BM
make bemHG_const_bm.mex.cpp
make bemHG_const_bm_sp.mex.cpp

% fmbem

make leafHp.mex.c
make shiftup.mex.c
make trans.mex.c
make recover.mex.c
make recover_bm.mex.c
make leafGq.mex.c
make shiftdown.mex.c

movefile(sprintf('*.%s', mexext), '../');

% make ibemB_const