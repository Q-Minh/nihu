function M = ibemRHS(mesh)

nNodes = size(mesh.Nodes,1);
elements = drop_IDs(mesh);
M = sparse(nNodes, nNodes);
[xi4, w4] = gaussquad2(5, 4);
[N4, dN4] = shapefun(xi4, 24);
[xi3, w3] = gaussquad2(5, 3);
[N3, dN3] = shapefun(xi3, 23);
nE = size(elements,1);
for e = 1 : nE
    d = mod(elements(e,2),10);
    elem = elements(e, 4+(1:d));
    nodes = mesh.Nodes(elem,2:4);
    switch d
        case 3
            dN = dN3;
            N = N3;
            w = w3;
        case 4
            dN = dN4;
            N = N4;
            w = w4;
    end
    jac = cross(dN(:,:,1) * nodes, dN(:,:,2) * nodes, 2);
    M(elem,elem) = M(elem,elem) + N' * diag(w.*sqrt(dot(jac,jac,2))) * N;
end
end
