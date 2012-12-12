function phat = solve_direct(K, M, q, Freq)

%% Solution in the frequency domain
nFreq = length(Freq);
phat = zeros(length(q), nFreq);
wb = waitbar(0, 'Direct solution');
for iFreq = 1 : nFreq
    freq = Freq(iFreq);
    om = 2*pi*freq;
    phat(:,iFreq) = (K - om^2*M) \ (-i*om*q);
    waitbar(iFreq/nFreq, wb);
end
close(wb);