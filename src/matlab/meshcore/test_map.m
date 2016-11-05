clear;

m = Mesh();
ps = PointSet(DefCoordinateSystem.C3);
ps.addPoints([0 0 0; 1 0 0; 1 1 0]);
%ps.addPoints([0 0 0; 1 0 0]);
es = ElementSet(ps, ShapeSet.LinearLine, uint32([1, 2]), uint32([1, 2; 2, 3]));

coords = es.map(es.ShapeSet.Domain.CornerNodes);