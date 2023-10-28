classdef DefCoordinateSystem < CoordinateSystem
    enumeration
        C1 (Space.D1, [0], zeros(0,1))
        C2 (Space.D2, [0 0], [1 0])
        C3 (Space.D3, [0 0 0], [1 0 0; 0 1 0])
    end
end
