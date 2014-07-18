function P = bempow(model, p, v)

if length(p) == size(model.Elements,1)
    [gcoord, gnorm, w] = geo2gauss(model, 1);
    I = .5 * real(p .* conj(v));
    P = I .* repmat(w, 1, size(I,2));
end