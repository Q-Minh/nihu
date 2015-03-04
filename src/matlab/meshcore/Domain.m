classdef Domain
    properties
        Space
        CornerNodes
    end
    
    methods
        function obj = Domain(space, cornernodes)
            obj.Space = space;
            obj.CornerNodes = cornernodes;
        end
    end
    
    enumeration
        Line(Space.D1, [-1; +1]);
        Tria(Space.D2, [0 0; 1 0; 0 1]);
        Quad(Space.D2, [-1 -1; 1 -1; 1 1; -1 1]);
    end
end
