function phat = solve_modal(Phi, Om, q, Freq)

%% Solution in the frequency domain
nFreq = length(Freq);
nModes = length(Om);
Lambda = diag(Om.^2);
I = eye(nModes);
alpha = zeros(nModes, nFreq);
wb = waitbar(0, 'Modal solution');
PhiTq = Phi.'*q;
for iFreq = 1 : nFreq
    freq = Freq(iFreq);
    om = 2*pi*freq;
    alpha(:,iFreq) = (Lambda - om^2*I) \ (-i*om*PhiTq);
    waitbar(iFreq/nFreq, wb);
end
close(wb);
phat = Phi * alpha;