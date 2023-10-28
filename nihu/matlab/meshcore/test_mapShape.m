clear;
D = Domain.Line;
[cs1, conn1] = D.divide(10);

[cs2, conn2] = mapShape(cs1, conn1, ShapeSet.LinearLine, ShapeSet.QuadraticLine);