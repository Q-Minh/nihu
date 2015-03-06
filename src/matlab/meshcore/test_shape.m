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
    disp(s.char());
    
    N = s.eval(s.Nodes);
    disp(N);
end
