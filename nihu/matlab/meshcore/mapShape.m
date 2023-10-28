function [newcs, newconn] = mapShape(cs, conn, oldShape, newShape)

N = oldShape.eval(newShape.Nodes);

for d = 1 : oldShape.Domain.Space.Dimension
    x = cs.coordinates(:,d);
    newxi(:,:,d) = N * x(conn');
end

N2 = size(newShape.Nodes,1);
nE = size(conn,1);

newconn = reshape(1 : nE * N2, N2, nE)';
newxi = reshape(newxi, [], oldShape.Domain.Space.Dimension);
newcs = CoordinateSet(newxi);
[~, ic] = newcs.mergeCoincident(1e-5);
newconn = ic(newconn);

end
