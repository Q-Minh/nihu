function y = bb_matvec_regular(x, T, P2P, P2M, M2M, M2L, L2L, L2P, alpha)
y = bb_far_transfer(T, x, P2M, M2M, M2L, L2L, L2P, alpha);
y = y + P2P * x;
end