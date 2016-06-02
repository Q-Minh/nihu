clear;
clc;
% clear classes;

%%
ps1 = PointSet(1, DefCoordinateSystem.C3);
ps1.addPoints(rand(100,3));

ps2 = PointSet(2, DefCoordinateSystem.C3);
ps2.addPoints(rand(100, 3), uint32(11:110));

newIds = ps1.merge(ps2);

ps1.addPoints([0 0 0], uint32(1000));
newIds = ps1.addPoints(rand(100,3));
p = ps1.getCoordinatesById(1:8);

e = [
    ShapeSet.LinearTria.Id 1 2 3
    ];

es = ElementSet(1, ps1);
es.addElem(e, 5);
es.addElem([ShapeSet.LinearLine.Id 25 2000]);
