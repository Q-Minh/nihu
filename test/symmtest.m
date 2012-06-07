clear;
clc;
close all;

kvec = 10 : 5 : 80;
nK = length(kvec);

t = zeros(nK, 4);
N = zeros(nK, 2);
for iK = 1 : nK
    t(iK, :) = symmetric_fuso(kvec(iK));
    N(iK, :) = symmetric_fuso2(kvec(iK));
end

figure;
plot(kvec, N);

figure;
plot(N(:,1), t(1:size(N,1),:));

figure;
tmatch = N(:,1).^2;
tmatch = t(end,4) ./ tmatch(end) * tmatch;
t2match = N(:,1).*log(N(:,1)).^2;
t2match = t(end,4) ./ t2match(end) * t2match;
plot(N(:,1),  t(:,4), N(:,1), tmatch, N(:,1), t2match);

save data/times2 t