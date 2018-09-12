function G = bem_matrix(mesh, kernel, singular)

nodes = mesh.Nodes(:,2:4);

[x, nx] = centnorm(mesh);

[xi, wns] = gaussquad2(40, 4);
[Nns, dNns] = ShapeSet.LinearQuad.eval(xi);

[xi, wr] = gaussquad2(6, 3);
[Nr, dNr] = ShapeSet.LinearTria.eval(xi);

nDof = size(mesh.Elements,1);

G = zeros(nDof, nDof);

for e = 1 : nDof
    
    progbar(1, nDof, e);
    
    xelem = mesh.Elements(e, 5:7);
    
    for f = 1 : nDof
        yelem = mesh.Elements(f, 5:7);
        
        if (e == f) % face match singular
            y = nodes(yelem,:);
            
            G(e,e) = singular(x(e,:), nx(e,:), y);
        elseif (any(ismember(yelem, xelem))) % nearly singular
            order = 40;
            
            
            y = Nns * nodes(yelem([1 2 3 3]),:);
            xgxi = dNns(:,:,1) * nodes(yelem([1 2 3 3]),:);
            xgeta = dNns(:,:,2) * nodes(yelem([1 2 3 3]),:);
            jvec = cross(xgxi, xgeta, 2);
            jac = sqrt(dot(jvec, jvec, 2));
            ny = bsxfun(@times, jvec, 1./jac);
            
            g = kernel(x(e,:), nx(e,:), y, ny);
            
            G(e,f) = wns.' * diag(jac) * g;
        else % regular integral
            order = 8;
            
            
            y = Nr * nodes(yelem,:);
            xgxi = dNr(:,:,1) * nodes(yelem,:);
            xgeta = dNr(:,:,2) * nodes(yelem,:);
            jvec = cross(xgxi, xgeta, 2);
            jac = sqrt(dot(jvec, jvec, 2));
            ny = bsxfun(@times, jvec, 1./jac);
            
            g = kernel(x(e,:), nx(e,:), y, ny);
            
            G(e,f) = wr.' * diag(jac) * g;
        end % of regular case
    end % of loop over cols
end % of loop over rows

end % of function