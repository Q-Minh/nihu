surf = create_sphere_boundary(1, 6);	%// radiating surface
field = create_sphere_boundary(1.2, 6);	%// field point mesh

[surf_nodes, surf_elem] = extract_Boonen_mesh(surf);
[field_nodes, field_elem] = extract_Boonen_mesh(field);

k = 5;	%// wave number 

%// call C++ program
[L, M, Mt, N, Lc, Mc, Mtc, Nc, t1, t2] =...
    couple_radiation(surf_nodes, surf_elem, field_nodes, field_elem, k);
    
%// compare result matrices and computation times
fprintf(1, 'L log10error:\t%.3f\n', log10(norm(Lc-L)/norm(Lc)));
fprintf(1, 'M log10error:\t%.3f\n', log10(norm(Mc-M)/norm(Mc)));
fprintf(1, 'Mt log10error:\t%.3f\n', log10(norm(Mtc-Mt)/norm(Mtc)));
fprintf(1, 'N log10error:\t%.3f\n', log10(norm(Nc-N)/norm(Nc)));
fprintf(1, 'time gain:\t%g\n', t2/t1);
