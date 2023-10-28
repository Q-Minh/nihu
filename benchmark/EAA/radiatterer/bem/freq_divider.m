c = 340;

flim = nan(5,2);
flim(1,1) = 300;

for d = 1 : 5
    Le = .1 / d;
    radi = create_radiatterer(Le);
    kmax = min(mesh_kmax(radi, 4))
    freq_max = kmax * c / 2/ pi;
    
    if d > 1
        flim(d,1) = round(flim(d-1,2) * .8);
    end
    flim(d,2) = round(freq_max * 1.17);
    
end

for d = 1 : 5
    Le = .1 / d;
    fprintf(1, 'inner_cycle_fmm %d %d %d %03d %g\n', ...
        flim(d,1)*2, 1, flim(d,2)*2, floor(100/d), 2*Le)
end
