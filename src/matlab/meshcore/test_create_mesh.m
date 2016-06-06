clear;

mesh = Mesh.createEmpty();
mesh.addElemSet(ElemSet.createLine([0, 0], [1, 2], 100, 'quadratic'));

pS = PointSet(DefCoordinateSystem.C3);
pS.addCoordinates( . )
eS = ElementSet(pS, 1, [ShapeSet.QuadraticLine.Id 1 2]);
eS.divide(100);