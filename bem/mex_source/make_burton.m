clc
clear mex

make bemHG_const_bm.mex.c

movefile(sprintf('*.%s', mexext), '../');

% make ibemB_const