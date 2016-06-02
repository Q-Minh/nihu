classdef DefCoordinateSystem < CoordinateSystem
    enumeration
        C2 (Space.D2, [0 0], [1 0])
        C3 (Space.D3, [0 0 0], [1 0 0; 0 1 0])
    end
end
