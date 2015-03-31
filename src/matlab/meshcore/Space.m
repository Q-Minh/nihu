classdef Space
    properties
        Dimension;
    end
    
    methods
        function obj = Space(d)
            obj.Dimension = d;
        end
    end
    
    enumeration
        D0(0)
        D1(1)
        D2(2)
        D3(3)
    end
end
