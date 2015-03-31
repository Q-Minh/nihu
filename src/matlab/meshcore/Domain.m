classdef Domain < handle
    properties (SetAccess = immutable)
        Space
        CornerNodes
        Volume
        Center
    end
    
    methods
        function obj = Domain(space, cornernodes, volume, center)
            obj.Space = space;
            obj.CornerNodes = cornernodes;
            obj.Volume = volume;
            obj.Center = center;
        end
    end
    
    enumeration
        Point(Space.D0, zeros(1,0), 0, zeros(1,0));
        Line(Space.D1, [-1; +1], 2, 0);
        Tria(Space.D2, [0 0; 1 0; 0 1], .5, [1 1]/3);
        Quad(Space.D2, [-1 -1; 1 -1; 1 1; -1 1], 4, [0 0]);
        Tetra(Space.D3, [0 0 0; 1 0 0; 0 1 0; 0 0 1], 1/6, [1 1 1]/4);
        Penta(Space.D3, [0 0 -1; 1 0 -1; 0 1 -1; 0 0 1; 1 0 1; 0 1 1], 1, [1 1 0]/3);
        Hexa(Space.D3, [
            -1 -1 -1; 1 -1 -1; 1 1 -1; -1 1 -1
            -1 -1 +1; 1 -1 +1; 1 1 +1; -1 1 +1
            ], 8, [0 0 0]);
    end
end
