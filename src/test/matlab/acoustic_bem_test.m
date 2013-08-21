clear;

surface = create_sphere_boundary(1, 5);
surface = quad2tria(surface);
field = create_sphere_boundary(3, 2);

x0 = [0 0 0];

[cs, ns] = centnorm(surface);
[cf, nf] = centnorm(field);
[nods, els] = extract_Boonen_mesh(surface);
[nodf, elf] = extract_Boonen_mesh(field);


kvec = pi-2e-2 : 1e-3 : pi+7e-2;
for iK = 1 : length(kvec)
    k = kvec(iK);
    disp(k);
    
    [ps_anal(:,iK), qs_anal(:,iK)] = incident('point', x0, cs, ns, k);
    [pf_anal(:,iK), qf_anal(:,iK)] = incident('point', x0, cf, nf, k);
    
    tic;
    [Gs, Hs, Hts, Ks, Gf, Hf, Htf, Kf, dur_sep] = acoustic_bem(nods, els, nodf, elf, k);
    toc;
    
    ps(:,iK) = Hs \ (Gs * qs_anal(:,iK));
    err_ps(:,iK) = log10(abs(ps(:,iK)./ps_anal(:,iK)-1));
    
    ps2(:,iK) = Ks \ (Hts * qs_anal(:,iK));
    err_ps2(:,iK) = log10(abs(ps2(:,iK)./ps_anal(:,iK)-1));
    
    alph = 1i / k;
    ps_bm(:,iK) = (Hs + alph * Ks) \ ((Gs + alph * Hts) * qs_anal(:,iK));
    err_ps_bm(:,iK) = log10(abs(ps_bm(:,iK)./ps_anal(:,iK)-1));
end

clear mex
