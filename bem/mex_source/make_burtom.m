clear mex

make bemHG_lin_bm.mex.c

movefile(sprintf('*.%s', mexext), '../');

% make ibemB_const