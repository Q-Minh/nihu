function M = ibemRHS(mesh)

nNodes = size(mesh.Nodes,1);
elements = drop_IDs(mesh);
M = sparse(nNodes, nNodes);

[xi3, w3] = gaussquad2(5, 3);
[N3, dN3] = shapefun(xi3, 23);

[xi4, w4] = gaussquad2(5, 4);
[N4, dN4] = shapefun(xi4, 24);

[xi42, w42] = gaussquad2(7, 4);
[N42, dN42] = shapefun(xi42, 242);

nE = size(elements,1);
for e = 1 : nE
    switch elements(e,2)
        case {23 231}
            elem = elements(e, 4+(1:3));
            nodes = mesh.Nodes(elem,2:4);
            dN = dN3;
            N = N3;
            w = w3;
        case {24 241}
            elem = elements(e, 4+(1:4));
            nodes = mesh.Nodes(elem,2:4);
            dN = dN4;
            N = N4;
            w = w4;
        case 242
            elem = elements(e, 4+(1:9));
            nodes = mesh.Nodes(elem,2:4);
            dN = dN42;
            N = N42;
            w = w42;
    end
    jac = cross(dN(:,:,1) * nodes, dN(:,:,2) * nodes, 2);
    M(elem,elem) = M(elem,elem) + N' * diag(w.*sqrt(dot(jac,jac,2))) * N;
end
end
