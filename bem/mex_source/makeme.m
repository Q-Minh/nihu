clear mex

make bemHG_const
make bemHG_lin
make ibemB_lin
make ibemD_const
make ibemD_lin

make bemHG_const_sp

make leafHp
make shiftup
make trans
make recover
make leafGq
make shiftdown

movefile(sprintf('*.%s', mexext), '../');

% make ibemB_const