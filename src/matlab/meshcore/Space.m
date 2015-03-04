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
        D1(1)
        D2(2)
        D3(3)
    end
end
