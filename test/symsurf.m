clear;
clc;

kvec = 10 : 5 : 40;
nK = length(kvec);
for iK = 1 : length(kvec)
    Time(iK,:,:) = symsurffun(kvec(iK));
end