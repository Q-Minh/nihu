clear;
m = Mesh();
ps = PointSet(DefCoordinateSystem.C3);
ps.addPoints([0 0 0; 1 0 0]);
ps.addPoints([0 0 0; 1 0 0]);
es = ElementSet(ps, uint32(1), [ShapeSet.LinearLine.Id, 1, 2]);
m.addElementSet(es);