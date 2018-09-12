function G = bem_matrix(mesh, kernel)

nodes = mesh.Nodes(:,2:4);

nDof = size(mesh.Elements,1);

G = zeros(nDof, nDof);

for e = 1 : nDof
    
    xelem = mesh.Elements(e, 5:7);
    
    for f = 1 : nDof
        yelem = mesh.Elements(f, 5:7);
        
        if (e == f) % face match singular
        elseif (any(ismember(yelem, xelem))) % nearly singular
            order = 40;
            
            [xi, w] = gaussquad2(order, 4);
            
            [N, dN] = ShapeSet.LinearQuad.eval(xi);
            y = N * nodes(elem,:);
            xgxi = dN(:,:,1) * nodes(yelem([1 2 3 3]),:);
            xgeta = dN(:,:,2) * nodes(yelem([1 2 3 3]),:);
            jvec = cross(xgxi, xgeta, 2);
            jac = sqrt(dot(jvec, jvec, 2));
            ny = bsxfun(@times, jvec, 1./jac);
            
            g = kernel(x, nx, y, ny);
            
            G(e,f) = w.' * diag(jac) * g;
        else % regular integral
            order = 8;
            
            [xi, w] = gaussquad2(order, 3);
            
            [N, dN] = ShapeSet.LinearTria.eval(xi);
            y = N * nodes(elem,:);
            xgxi = dN(:,:,1) * nodes(yelem,:);
            xgeta = dN(:,:,2) * nodes(yelem,:);
            jvec = cross(xgxi, xgeta, 2);
            jac = sqrt(dot(jvec, jvec, 2));
            ny = bsxfun(@times, jvec, 1./jac);
            
            g = kernel(x, nx, y, ny);
            
            G(e,f) = w.' * diag(jac) * g;
        end % of regular case
    end % of loop over cols
end % of loop over rows

end % of function