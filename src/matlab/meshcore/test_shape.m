clear;
clc;

shapes = [
    ShapeSet.LinearLine
    ShapeSet.QuadraticLine
    
    ShapeSet.LinearTria
    ShapeSet.QuadraticTria
    
    ShapeSet.LinearQuad
    ShapeSet.QuadraticQuad
    
    ShapeSet.LinearTetra
    ShapeSet.LinearPenta
    ShapeSet.LinearHexa
    ];

for i = 1 : length(shapes)
    s = shapes(i);
    fprintf('Testing Domain %s...\n', s.char());
    
    %% check if nodal shape functions are identity matrix
    N = s.eval(s.Nodes);
    err = norm(N-eye(size(N)), 'fro')/norm(eye(size(N)), 'fro');
    if (err > eps)
        error('ShapeSet identity error: %s', s.char());
    end
    
    %%
    disp('Tests OK');
end

