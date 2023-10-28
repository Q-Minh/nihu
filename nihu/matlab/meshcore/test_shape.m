clear;
clc;

shapes = [
    ShapeSet.ConstantLine
    ShapeSet.LinearLine
    ShapeSet.QuadraticLine
    
    ShapeSet.ConstantTria
    ShapeSet.LinearTria
    ShapeSet.QuadraticTria
    
    ShapeSet.ConstantQuad
    ShapeSet.LinearQuad
    ShapeSet.QuadraticQuad
    
    ShapeSet.ConstantTetra
    ShapeSet.LinearTetra
    ShapeSet.ConstantPenta
    ShapeSet.LinearPenta
    ShapeSet.ConstantHexa
    ShapeSet.LinearHexa
    ];

for i = 1 : length(shapes)
    s = shapes(i);
    fprintf('Testing ShapeSet %s...\n', s.char());
    
    %% check if nodal shape functions are identity matrix
    N = s.eval(s.Nodes);
    err = norm(N-eye(size(N)), 'fro')/norm(eye(size(N)), 'fro');
    if (err > eps)
        error('ShapeSet identity error: %s', s.char());
    end
    
    %%
    disp('Tests OK');
end

