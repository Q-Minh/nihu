mex M2L.cpp -ID:/Peti/research/toolbox/eigen -output compute_M2L

%%
multi = ones(64,12);
s = [1 2];
r = [1 2];
ridx = [1 2];
local = compute_M2L(int32(s)-1, int32(r)-1, int32(ridx)-1, M2L, multi);