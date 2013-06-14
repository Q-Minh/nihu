clear;
for i = 1 : 13
    disp(i);
    
    [XI3{i}, W3{i}] = edgematch_tria(i);
    sum(W3{i})
    I3(i) = integrate(XI3{i}, W3{i});
    n3(i) = size(XI3{i},1);

    [XI34{i}, W34{i}] = edgematch_mixed(i);
    sum(W34{i})
    I34(i) = integrate(XI34{i}, W34{i});
    n34(i) = size(XI34{i},1);

    [XI4{i}, W4{i}] = edgematch_quad(i);
    I4(i) = integrate(XI4{i}, W4{i});
    n4(i) = size(XI4{i},1);

    [XIG{i}, WG{i}] = edgematch_quad_gl(i);
    IG(i) = integrate(XIG{i}, WG{i});
    sum(WG{i})
    nG(i) = size(XIG{i},1);

end

%%
figure;
hold on;
semilogx(n3, log10(abs(I3/I3(end)-1)), 'b*-');
semilogx(n4, log10(abs(I4/I4(end)-1)), 'r*-');
semilogx(n34, log10(abs(I34/I34(end)-1)), 'mo-');
semilogx(nG, log10(abs(IG/I4(end)-1)), 'ko-');
set(gca, 'xscale', 'log');
