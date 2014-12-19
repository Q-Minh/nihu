function P2P = bb_P2P_regular(leaflevel, xRec, xSou, sp_kernel)
%BB_P2P_REGULAR compute P2P sparse matrix of a bb FMM

[i, j] = nfij(leaflevel.nearfield, leaflevel.nodsou, leaflevel.nodrec);
z = sp_kernel(xRec(i,:), xSou(j,:));
P2P = sparse(i, j, z);

end
