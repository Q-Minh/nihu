clear;
cs = CoordinateSet;

k = 0;
while true
    k = k+1;
    idx = cs.addCoordinates(rand(1000,3), 'merge', 1e-3);
    if max(idx) ~= k * 1000
        break;
    end
end
