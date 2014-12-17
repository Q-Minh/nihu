mex M2L.cpp -I/D/research/toolbox/eigen -output compute_M2L
mex M2M.cpp -I/D/research/toolbox/eigen -output compute_M2M
mex L2L.cpp -I/D/research/toolbox/eigen -output compute_L2L

movefile *.mexa* ../