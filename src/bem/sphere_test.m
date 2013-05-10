clear;

for is_quadratic = [false true]
    for N = 4 : 15
        for tria_num = [0 1 2] % 0 quad 1 half 2 full tria
            for is_isoparam = [false true]
                
                disp([is_quadratic N is_isoparam tria_num]);
                
                try
                    clear G H M
                    
                    R = 1;                                  % radius
                    mesh = create_sphere_boundary(1, N);    % mesh
                    k = min(mesh_kmax(mesh)) / 2;           % wave number
                    if tria_num == 2
                        mesh = quad2tria(mesh);
                    elseif tria_num == 1
                        eps = 1e-3;
                        mup = mesh_section(mesh, [-Inf, Inf; -eps Inf; -Inf Inf]');
                        mlo = mesh_section(mesh, [-Inf, Inf; -Inf eps; -Inf Inf]');
                        mlo = quad2tria(mlo);
                        mesh = join_meshes(mup, mlo);
                        mesh = merge_coincident_nodes(mesh);
                        mesh = drop_unused_nodes(mesh);
                    end
                    [~, ~, w] = geo2gauss(mesh, 0);
                    
                    if is_quadratic
                        mesh = quadratise(mesh);
                        r = sqrt(dot(mesh.Nodes(:,2:4), mesh.Nodes(:,2:4), 2));
                        mesh.Nodes(:,2:4) = bsxfun(@times, mesh.Nodes(:,2:4), R./r);
                    end
                    
                    [nodes, elements] = extract_Boonen_mesh(mesh);
                    
                    if is_isoparam
                        tic;
                        [G, H, neval] = Boonen13_gal_iso(nodes, elements, k);
                        tBooni = toc;
                        M = ibemRHS(mesh);
                    else
                        tic;
                        [G, H, neval] = Boonen13_gal_const(nodes, elements, k);
                        tBooni = toc;
                        M = diag(w);
                    end
                    
                    p = (H - .5*M) \ sum(G,2);
                    p_anal = -R/(1 + 1i*k*R);
                    
                    save(fullfile('data', sprintf('res_N%u_q%u_f%u_t%u', N, is_quadratic, is_isoparam, tria_num)), ...
                        'p', 'p_anal', 'tBooni', 'neval');
                    
                end
            end
        end
    end
end
